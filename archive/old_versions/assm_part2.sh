#!/bin/bash
set -euo pipefail

# --- Usage information ---
show_help() {
    echo "Usage: $0 --input <input_file> --assembler <flye|canu> --genomeSize <size> [options]"
    echo
    echo "Required:"
    echo "  --input <file>        Input metadata file"
    echo "  --genomeSize <size>   Expected genome size (e.g., 40m or 2.5g. Only vital if you are assembling with canu! Otherwise, you can put anything here.)"
    echo
    echo "Optional:"
    echo "  --threads <N>         Number of threads (auto-detected if not provided)"
    echo "  --maxMemory <GB>      Max memory (auto-detected if not provided)"
    exit 1
}

# --- Parse arguments ---
INPUT_FILE=""
GENOME_SIZE=""
MAX_MEMORY=""
THREADS="${SLURM_CPUS_ON_NODE:-$(nproc)}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            [[ $# -ge 2 ]] || { echo "Error: --input requires a file argument"; exit 1; }
            INPUT_FILE="$2"; shift 2 ;;
        --genomeSize)
            [[ $# -ge 2 ]] || { echo "Error: --genomeSize requires a value"; exit 1; }
            GENOME_SIZE="$2"; shift 2 ;;
        --threads)
            [[ $# -ge 2 ]] || { echo "Error: --threads requires a value"; exit 1; }
            THREADS="$2"; shift 2 ;;
        --maxMemory)
            [[ $# -ge 2 ]] || { echo "Error: --maxMemory requires a value"; exit 1; }
            MAX_MEMORY="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Unknown option: $1"; show_help ;;
    esac
done

# --- Validate arguments ---
if [[ -z "$INPUT_FILE" || -z "$GENOME_SIZE" ]]; then
    echo "Error: --input and --genomeSize are required."
    show_help
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# --- Extract OUTPUT_DIR early (set -e safe) ---
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Could not extract Desired Output Location from $INPUT_FILE"
    exit 1
fi

# --- Auto-detect memory if not provided ---
if [[ -z "${MAX_MEMORY:-}" ]]; then
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        MAX_MEMORY=$((SLURM_MEM_PER_NODE / 1024))
    elif command -v free >/dev/null 2>&1; then
        MAX_MEMORY=$(free -g | awk '/^Mem:/{print int($2*0.8)}')
    else
        MAX_MEMORY=128
    fi
fi

echo "Detected configuration:"
echo "  Input File:  $INPUT_FILE"
echo "  Threads:     $THREADS"
echo "  Memory:      ${MAX_MEMORY} GB"
echo "  Genome Size: $GENOME_SIZE"
echo

# --- Helper: check that a file exists and is not empty ---
check_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -s "$file" ]]; then
        echo "Error: $desc ($file) not found or empty."
        exit 1
    fi
}

# --- Define working directories ---
LOGDIR="$OUTPUT_DIR/logs"
ASM_DIR="$OUTPUT_DIR/assembled_genomes"
mkdir -vp "$OUTPUT_DIR" "$LOGDIR" "$ASM_DIR"

# --- Sanitized DES collisions ---
declare -A DES_MAP
declare -i DES_COLLISION=0

while IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER BUSCO_LIN; do
    [[ -z "${LR:-}" ]] && continue
    [[ "$LR" =~ ^# ]] && continue
    [[ -z "${DES:-}" ]] && continue

    BASE="${DES//[^a-zA-Z0-9._-]/_}"

    if [[ -n "${DES_MAP[$BASE]:-}" ]]; then
        if [[ "${DES_MAP[$BASE]}" != *"|$DES|"* ]]; then
            DES_MAP[$BASE]+="$DES|"
        fi
        DES_COLLISION=1
    else
        DES_MAP[$BASE]="|$DES|"
    fi
done < <(sed 's/\r$//' "$INPUT_FILE")

if (( DES_COLLISION )); then
    echo "Error: DES collisions detected after sanitization:"
    for base in "${!DES_MAP[@]}"; do
        IFS='|' read -ra DES_LIST <<< "${DES_MAP[$base]}"
        if (( ${#DES_LIST[@]} > 2 )); then
            echo "  '$base' <= ${DES_LIST[*]:1:-1}"
        fi
    done
    exit 1
fi

unset DES_MAP DES_COLLISION

# --- Process each entry in the metadata file ---
while IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER BUSCO_LIN; do
    [[ -z "${LR:-}" ]] && continue
    [[ "$LR" =~ ^# ]] && continue

    if [[ -z "$LR" || -z "$R1" || -z "$R2" || -z "$RNA1" || -z "$RNA2" || -z "$DES" || -z "$SPECIES" || -z "$ASSEMBLER" || -z "$BUSCO_LIN" ]]; then
        echo "Warning: Skipping malformed line: $LR $R1 $R2 $RNA1 $RNA2 $DES $SPECIES $ASSEMBLER $BUSCO_LIN"
        continue
    fi

    ASSEMBLER=$(echo "$ASSEMBLER" | tr '[:upper:]' '[:lower:]')

    if [[ "$ASSEMBLER" != "flye" && "$ASSEMBLER" != "canu" && "$ASSEMBLER" != "spades" ]]; then
        echo "Error: Assembler must be 'flye', 'canu', or 'spades'."
        exit 1
    fi

    if [[ "$ASSEMBLER" == "spades" ]]; then
        command -v AAFTF >/dev/null 2>&1 || { echo >&2 "Error: AAFTF is not installed or not in PATH."; exit 1; }
        INPUT_DIR=$(sed -n 's/^#Location of Short Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        TYPE="SHORT"
    else
        command -v "$ASSEMBLER" >/dev/null 2>&1 || { echo >&2 "Error: $ASSEMBLER not found in PATH."; exit 1; }
        INPUT_DIR=$(sed -n 's/^#Location of Long Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
        TYPE="LONG"
    fi

    if [[ "$INPUT_DIR" == "NA" || -z "$INPUT_DIR" ]]; then
        echo "Error: Required input directory not specified for $ASSEMBLER."
        exit 1
    fi

    if [[ ! -d "$INPUT_DIR" ]]; then
        echo "Error: Input reads directory '$INPUT_DIR' does not exist."
        exit 1
    fi

    BASE="${DES//[^a-zA-Z0-9._-]/_}"

    if [[ "$TYPE" == "LONG" && "$LR" == "NA" ]]; then
        echo "No long-read file listed for $DES. Skipping."
        continue
    fi

    if [[ "$TYPE" == "SHORT" && ( "$R1" == "NA" || "$R2" == "NA" ) ]]; then
        echo "No short-read files listed for $DES. Skipping."
        continue
    fi

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part2.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part2.$BASE.$timestamp.log"

    if [[ "$TYPE" == "LONG" ]]; then
        INFILE="$INPUT_DIR/$LR"
        check_file "$INFILE" "Input long reads"

        if [[ "$ASSEMBLER" == "flye" ]]; then
            ASMFILE="$ASM_DIR/$BASE.flyeout/assembly.fasta"
        else
            ASMFILE="$ASM_DIR/$BASE.canuout/$BASE.contigs.fasta"
        fi

        {
            echo "--- Assembly for $DES using $ASSEMBLER ---"
            echo "Started at: $(date)"

            if [[ ! -s "$ASMFILE" ]]; then
                if [[ "$ASSEMBLER" == "flye" ]]; then
                    flye --nano-raw "$INFILE" --threads "$THREADS" -o "$ASM_DIR/$BASE.flyeout"
                else
                    canu -p "$BASE" -d "$ASM_DIR/$BASE.canuout" -nanopore "$INFILE" \
                        useGrid=false maxThreads="$THREADS" maxMemory="$MAX_MEMORY" \
                        stopOnLowCoverage=0 corOutCoverage=100 genomeSize="$GENOME_SIZE"
                fi
                check_file "$ASMFILE" "$ASSEMBLER assembly"
            else
                echo "Assembly already exists: $ASMFILE"
            fi

            echo "Completed at: $(date)"
        } > "$log_temp" 2>&1

    else
        ASMFILE="$ASM_DIR/$BASE.spades.fasta"
        {
            echo "--- Assembly for $DES using spades ---"
            echo "Started at: $(date)"

            check_file "$INPUT_DIR/$R1" "Left short reads"
            check_file "$INPUT_DIR/$R2" "Right short reads"

            if [[ ! -s "$ASMFILE" ]]; then
                AAFTF assemble --isolate -c "$THREADS" \
                    --left "$INPUT_DIR/$R1" --right "$INPUT_DIR/$R2" \
                    --memory "$MAX_MEMORY" -o "$ASMFILE" -w "$ASM_DIR/$BASE.spadesout"
                check_file "$ASMFILE" "SPAdes assembly"
                rm -rf "$ASM_DIR/$BASE.spadesout/K??" "$ASM_DIR/$BASE.spadesout/tmp"
            else
                echo "Assembly already exists: $ASMFILE"
            fi

            echo "Completed at: $(date)"
        } > "$log_temp" 2>&1
    fi

    sed -E "s/\x1B\[[0-9;]*[mK]//g" "$log_temp" > "$log_file"
    rm -f "$log_temp"

done < <(sed 's/\r$//' "$INPUT_FILE")

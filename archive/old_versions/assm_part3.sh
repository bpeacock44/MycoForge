#!/bin/bash
set -euo pipefail

# --- Usage information ---
show_help() {
    echo "Usage: $0 --input <input_file> --assembler <flye|canu> [--threads <N>] [--maxMemory <GB>]"
    echo
    echo "Required arguments:"
    echo "  --input <file>          Metadata input file"
    echo
    echo "Optional arguments:"
    echo "  --threads <N>           Number of threads (auto-detected if not provided)"
    echo "  --maxMemory <GB>        Max memory in GB (auto-detected from SLURM if not provided)"
    exit 1
}

# --- Parse arguments ---
INPUT_FILE=""
MAX_MEMORY=""
THREADS="${SLURM_CPUS_ON_NODE:-$(nproc)}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            [[ $# -ge 2 ]] || { echo "Error: --input requires a file argument"; exit 1; }
            INPUT_FILE="$2"; shift 2 ;;
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

# --- Validate required inputs ---
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: --input is required."
    show_help
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# --- Check dependencies ---
for cmd in AAFTF minimap2 medaka_consensus; do
    command -v "$cmd" >/dev/null 2>&1 || { echo >&2 "Error: $cmd not found in PATH."; exit 1; }
done

# --- Auto-detect memory if not provided ---
if [[ -z "${MAX_MEMORY:-}" ]]; then
    if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
        MAX_MEMORY=$((SLURM_MEM_PER_NODE / 1024))  # SLURM gives KB
    elif command -v free >/dev/null 2>&1; then
        MAX_MEMORY=$(free -g | awk '/^Mem:/{print int($2*0.8)}')
    else
        MAX_MEMORY=128  # Fallback
    fi
fi

echo "Detected configuration:"
echo "  Input file: $INPUT_FILE"
echo "  Threads:    $THREADS"
echo "  Max Memory: ${MAX_MEMORY} GB"
echo

# --- Extract directories from header ---
LR_DIR=$(sed -n 's/^#Location of Long Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
SR_DIR=$(sed -n 's/^#Location of Short Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$LR_DIR" || -z "$SR_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Could not extract required directories from $INPUT_FILE"
    exit 1
fi

# --- Validate input directory ---
[[ "$LR_DIR" != "NA" ]] && [[ ! -d "$LR_DIR" ]] && echo "Error: $LR_DIR does not exist." && exit 1
[[ "$SR_DIR" != "NA" ]] && [[ ! -d "$SR_DIR" ]] && echo "Error: $SR_DIR does not exist." && exit 1

# --- Helper to verify file existence ---
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
WORKDIR="$OUTPUT_DIR/processed_small_reads"
ASM_DIR="$OUTPUT_DIR/assembled_genomes"

mkdir -vp "$LOGDIR" "$WORKDIR" "$ASM_DIR"

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

# Example header in file:
# #Location of Long Reads=/path/to/long_read_data
# #Location of Short Reads=/path/to/short_read_data
# #Location of RNAseq Reads=/path/to/RNAseq_read_data
# #Desired Output Location=/path/to/output
# #LR R1 R2 RNA1 RNA2 Designation Species_Name Assembler Busco_Lineage

# --- Process each entry in the metadata file ---
# normalize CRLFs via sed and feed into the while loop
while IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER BUSCO_LIN; do
    # skip empty lines
    [[ -z "${LR:-}" ]] && continue
    [[ "$LR" =~ ^# ]] && continue

    # skip malformed lines
    if [[ -z "$LR" || -z "$R1" || -z "$R2" || -z "$RNA1" || -z "$RNA2" || -z "$DES" || -z "$SPECIES" || -z "$ASSEMBLER" || -z "$BUSCO_LIN" ]]; then
        echo "Warning: Skipping malformed line: $LR $R1 $R2 $RNA1 $RNA2 $DES $SPECIES $ASSEMBLER $BUSCO_LIN"
        continue
    fi

    # normalize to lowercase
    ASSEMBLER=$(echo "$ASSEMBLER" | tr '[:upper:]' '[:lower:]')
    
    # make sure assembler is valid
    if [[ "$ASSEMBLER" != "flye" && "$ASSEMBLER" != "canu" && "$ASSEMBLER" != "spades" ]]; then
        echo "Error: Assembler must be 'flye', 'canu', or 'spades'."
        exit 1
    fi

    LONG_FLAG=0
    [[ "$LR" != "NA" ]] && LONG_FLAG=1

    SHORT_FLAG=0
    [[ "$R1" != "NA" ]] && SHORT_FLAG=1

    # --- Define file paths #
    BASE=${DES//[^a-zA-Z0-9._-]/_}

    ONT=""
    (( LONG_FLAG )) && ONT="$LR_DIR/$LR"
    SR_R1="$SR_DIR/$R1"
    LEFT="$WORKDIR/${BASE}_filtered_1.fastq.gz"
    RIGHT="$WORKDIR/${BASE}_filtered_2.fastq.gz"

    POLISHDIR="$ASM_DIR/$BASE.$ASSEMBLER.polishout"
    mkdir -vp "$POLISHDIR"

    if [[ "$ASSEMBLER" == "flye" ]]; then
        ASMFILE="$ASM_DIR/$BASE.flyeout/assembly.fasta"
    elif [[ "$ASSEMBLER" == "canu" ]]; then
        ASMFILE="$ASM_DIR/$BASE.canuout/$BASE.contigs.fasta"
    else 
        ASMFILE="${ASM_DIR}/$BASE.spades.fasta"
    fi

    VECCLEAN="$POLISHDIR/$BASE.vecscreen.fasta"
    MEDAKA_OUT="$POLISHDIR/$BASE.medakaout"
    if (( LONG_FLAG )); then
        mkdir -p "$MEDAKA_OUT"
    fi
    MEDAKA_ASM="$MEDAKA_OUT/consensus.fasta"
    MEDAKA_MODEL="r1041_e82_400bps_hac_v5.0.0"
    ONT_PAF="$MEDAKA_OUT/$BASE.paf"
    POLCA="$POLISHDIR/$BASE.polca.fasta"
    SORTED="$POLISHDIR/$BASE.sorted.fasta"
    STATS="$POLISHDIR/$BASE.sorted.stats.txt"

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part3.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part3.$BASE.$timestamp.log"

    {
        echo "--- Post-assembly processing for $DES ($ASSEMBLER) ---"
        echo "Input file: $INPUT_FILE | Threads: $THREADS | Memory: $MAX_MEMORY GB"
        echo "Started at: $(date)"
        echo
        check_file "$ASMFILE" "Assembly"
        echo
        echo "--- Step 1: Vecscreen ---"
        if [[ ! -s "$VECCLEAN" ]]; then
            max_retries=3
            for i in $(seq 1 $max_retries); do
                echo "Running Vecscreen (attempt $i/$max_retries)..."
                if AAFTF vecscreen -i "$ASMFILE" -c "$THREADS" -o "$VECCLEAN"; then
                    echo "Vecscreen completed successfully."
                    break
                else
                    echo "Vecscreen failed (attempt $i/$max_retries)."
                    if [[ $i -lt $max_retries ]]; then
                        echo "Retrying in 30s..."
                        sleep 30
                    else
                        echo "Vecscreen failed after $max_retries attempts. Exiting."
                        exit 1
                    fi
                fi
            done
        else
            echo "$VECCLEAN exists; skipping vecscreen."
        fi
        check_file "$VECCLEAN" "Vecscreen output"

        echo
        echo "--- Step 2: Long-read polishing (Medaka) ---"
        if (( LONG_FLAG )); then
            if [[ ! -s "$MEDAKA_ASM" ]]; then
                minimap2 -t "$THREADS" -x map-ont "$VECCLEAN" "$ONT" > "$ONT_PAF"
                medaka_consensus -i "$ONT" -d "$VECCLEAN" -o "$MEDAKA_OUT" \
                                -t "$THREADS" -m "$MEDAKA_MODEL"
            else
                echo "Medaka output exists; skipping."
            fi
            check_file "$MEDAKA_ASM" "Medaka consensus output"
        else
            echo "No long-read data designated. Skipping long-read polishing (Medaka)"
        fi

        echo
        echo "--- Step 3: Short-read polishing (Polca) ---"
        if (( SHORT_FLAG )); then
            [[ $LONG_FLAG -eq 0 ]] && MEDAKA_ASM="$ASMFILE"
            if [[ ! -s "$POLCA" ]]; then
                AAFTF polish --method polca -i "$MEDAKA_ASM" -o "$POLCA" \
                            -c "$THREADS" --left "$LEFT" --right "$RIGHT" --memory "$MAX_MEMORY"
            else
                echo "$POLCA exists; skipping Polca."
            fi
            check_file "$POLCA" "Polca-polished assembly"
        else
            echo "No short-read data designated. Skipping short-read polishing (Polca)"
        fi

        echo
        echo "--- Step 4: Sorting assembly ---"
        # Decide what to sort:
        if [[ -s "$POLCA" ]]; then
            SORT_INPUT="$POLCA"
        elif [[ -s "$MEDAKA_ASM" ]]; then
            SORT_INPUT="$MEDAKA_ASM"
        else
            echo "No Polca or Medaka assembly available to sort for $DES. Skipping sort."
            SORT_INPUT=""
        fi
        
        if [[ -n "$SORT_INPUT" ]]; then
            if [[ ! -s "$SORTED" ]]; then
                echo "Sorting assembly: $SORT_INPUT"
                AAFTF sort -i "$SORT_INPUT" -o "$SORTED"
            else
                echo "$SORTED exists; skipping sort."
            fi
            check_file "$SORTED" "Sorted assembly"
        fi

        echo
        echo "--- Step 5: Assessing final assembly ---"
        if [[ ! -s "$STATS" ]]; then
            AAFTF assess -i "$SORTED" -r "$STATS"
        else
            echo "$STATS exists; skipping assess."
        fi

        echo
        echo "--- Finished post-assembly for $DES ---"
    } > "$log_temp" 2>&1

    sed -E "s/\x1B\[[0-9;]*[mK]//g" "$log_temp" > "$log_file"
    rm -f "$log_temp"

done < <(sed 's/\r$//' "$INPUT_FILE")

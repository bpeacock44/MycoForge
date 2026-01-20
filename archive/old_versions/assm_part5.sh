#!/bin/bash
set -euo pipefail


# --- Usage information ---
show_help() {
    echo "Usage: $0 --input <input_file> --type <protein|genome|both> [options]"
    echo
    echo "Optional arguments:"
    echo "  --threads <N>    Number of threads (auto-detects from SLURM or nproc)"
    exit 1
}

# --- Parse arguments ---
INPUT_FILE=""
TYPE=""
THREADS="${SLURM_CPUS_ON_NODE:-$(nproc)}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            [[ $# -ge 2 ]] || { echo "Error: --input requires a file argument"; exit 1; }
            INPUT_FILE="$2"; shift 2 ;;
        --type)
            [[ $# -ge 2 ]] || { echo "Error: --type requires a value"; exit 1; }
            TYPE="$2"; shift 2 ;;
        --threads)
            [[ $# -ge 2 ]] || { echo "Error: --threads requires a value"; exit 1; }
            THREADS="$2"; shift 2 ;;
        -h|--help) show_help ;;
        *) echo "Unknown option: $1"; show_help ;;
    esac
done

# --- Validate required inputs ---
if [[ -z "$INPUT_FILE" || -z "$TYPE" ]]; then
    echo "Error: Missing required arguments."
    show_help
fi

if [[ "$TYPE" != "protein" && "$TYPE" != "genome" && "$TYPE" != "both" ]]; then
    echo "Error: BUSCO type must be 'protein', 'genome', or 'both'."
    exit 1
fi

# --- Check dependencies ---
for cmd in busco; do
    command -v "$cmd" >/dev/null 2>&1 || { echo >&2 "Error: $cmd not found in PATH."; exit 1; }
done

echo "Detected configuration:"
echo "  Type:       $TYPE"
echo "  Threads:    $THREADS"
echo

# --- Extract directories from header ---
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Could not extract desired output directory from $INPUT_FILE"
    exit 1
fi

# --- Define working directories ---
LOGDIR="$OUTPUT_DIR/logs"
ANNO_DIR="$OUTPUT_DIR/funannotate_genomes"
BUSCO_DIR="$OUTPUT_DIR/busco_output"

mkdir -vp "$BUSCO_DIR" "$LOGDIR"

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

    # --- Define file paths ---
    BASE=${DES//[^a-zA-Z0-9._-]/_}

    PROT="$ANNO_DIR/$BASE.$ASSEMBLER/$BASE.funannotate04/annotate_results/${SPECIES}_${BASE}.proteins.fa"
    GEN="$ANNO_DIR/$BASE.$ASSEMBLER/$BASE.funannotate04/annotate_results/${SPECIES}_${BASE}.contigs.fsa"

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part5.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part5.$BASE.$timestamp.log"

    {
        echo "--- BUSCO for $BASE ($ASSEMBLER) ---"
        echo "Input file: $INPUT_FILE | Threads: $THREADS"
        echo "Busco lineage: $BUSCO_LIN"
        echo "Mode(s): $TYPE"
        echo

        # --- Run BUSCO on protein file ---
        if [[ "$TYPE" == "protein" || "$TYPE" == "both" ]]; then
            if [[ -s "$PROT" ]]; then
                echo "Running BUSCO (protein mode) on $PROT"
                busco -i "$PROT" \
                    --out_path "$BUSCO_DIR/$BASE.$ASSEMBLER.$BUSCO_LIN.protein.buscoout" \
                    -l "$BUSCO_LIN" \
                    -m protein \
                    -c "$THREADS"
            else
                echo "Warning: Protein file not found for $BASE: $PROT"
            fi
        fi

        # --- Run BUSCO on genome file ---
        if [[ "$TYPE" == "genome" || "$TYPE" == "both" ]]; then
            if [[ -s "$GEN" ]]; then
                echo "Running BUSCO (genome mode) on $GEN"
                busco -i "$GEN" \
                    --out_path "$BUSCO_DIR/$BASE.$ASSEMBLER.$BUSCO_LIN.genome.buscoout" \
                    -l "$BUSCO_LIN" \
                    -m genome \
                    -c "$THREADS"
            else
                echo "Warning: Genome file not found for $BASE: $GEN"
            fi
        fi

        echo "--- BUSCO completed for $BASE ---"
    } > "$log_temp" 2>&1

    # --- Strip ANSI color codes and finalize log ---
    sed -E "s/\x1B\[[0-9;]*[mK]//g" "$log_temp" > "$log_file"
    rm -f "$log_temp"

done < <(sed 's/\r$//' "$INPUT_FILE")
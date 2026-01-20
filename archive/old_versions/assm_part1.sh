#!/bin/bash
set -euo pipefail

# --- Usage information ---
show_help() {
    echo "Usage: $0 --input <input_file>"
    echo
    echo "Required:"
    echo "  --input <file>        Input metadata file"
    exit 1
}

# --- Parse arguments ---
INPUT_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            [[ $# -ge 2 ]] || { echo "Error: --input requires a file argument"; exit 1; }
            INPUT_FILE="$2"
            shift 2
            ;;
        -h|--help) show_help ;;
        *) echo "Unknown option: $1"; show_help ;;
    esac
done

# --- Check AAFTF dependency ---
command -v AAFTF >/dev/null 2>&1 || { echo >&2 "Error: AAFTF is not installed or not in PATH. Exiting."; exit 1; }

# --- Validate input ---
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: --input <file> is required."
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# --- Extract directories from header (set -e safe) ---
INPUT_DIR=$(sed -n 's/^#Location of Short Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Missing Short Read or Output directory in $INPUT_FILE"
    exit 1
fi

# --- Validate input directory ---
if [[ "$INPUT_DIR" == "NA" ]]; then
    echo "No short-read input directory designated. Skipping short-read preprocessing."
    exit 0
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Input reads directory '$INPUT_DIR' does not exist."
    exit 1
fi

# --- Sanitized DES collisions ---
declare -A DES_MAP
while IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER BUSCO_LIN; do
    [[ -z "${LR:-}" ]] && continue
    [[ "$LR" =~ ^# ]] && continue
    [[ -z "${DES:-}" ]] && continue

    BASE="${DES//[^a-zA-Z0-9._-]/_}"

    if [[ -n "${DES_MAP[$BASE]:-}" ]]; then
        echo "Error: Two or more designations are too similar when simplified. Please make them more unique."
        echo "  '$DES' and '${DES_MAP[$BASE]}' both map to '$BASE'"
        exit 1
    fi

    DES_MAP[$BASE]="$DES"
done < <(sed 's/\r$//' "$INPUT_FILE")

unset DES_MAP


# --- Helper: check that a file exists and is not empty ---
check_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -s "$file" ]]; then
        echo "Error: $desc ($file) not found or empty. Exiting."
        exit 1
    fi
}

# --- Define working directories ---
CPU=${SLURM_CPUS_ON_NODE:-1}
WORKDIR="$OUTPUT_DIR/processed_small_reads"
LOGDIR="$OUTPUT_DIR/logs"
mkdir -vp "$OUTPUT_DIR" "$WORKDIR" "$LOGDIR"

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

    if [[ "$R1" == "NA" || "$R2" == "NA" ]]; then
        timestamp=$(date +%Y%m%d_%H%M%S)
        log_file="$LOGDIR/part1.${DES}.${timestamp}.log"
        {
            echo "--- Skipping short-read preprocessing for $DES ---"
            echo "Reason: No short-read files listed in metadata (R1 or R2 = NA)"
            echo "Metadata file: $INPUT_FILE"
            echo "Time: $(date)"
        } | tee "$log_file"
        continue
    fi

    INFILE1="$INPUT_DIR/$R1"
    INFILE2="$INPUT_DIR/$R2"

    BASE="${DES//[^a-zA-Z0-9._-]/_}"

    LEFTTRIM="$WORKDIR/${BASE}_1P.fastq.gz"
    RIGHTTRIM="$WORKDIR/${BASE}_2P.fastq.gz"
    LEFT="$WORKDIR/${BASE}_filtered_1.fastq.gz"
    RIGHT="$WORKDIR/${BASE}_filtered_2.fastq.gz"

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp=$(mktemp "$LOGDIR/part1.${BASE}.XXXXXX.log.tmp")
    log_file="$LOGDIR/part1.${BASE}.${timestamp}.log"

    {
        echo "--- Preprocessing for $DES ---"
        echo "Metadata file: $INPUT_FILE"
        echo "Started at: $(date)"
        check_file "$INFILE1" "Input left reads"
        check_file "$INFILE2" "Input right reads"

        echo
        echo "--- Step 1: Trimming reads ---"
        if [[ ! -s "$LEFTTRIM" || ! -s "$RIGHTTRIM" ]]; then
            AAFTF trim --method bbduk --left "$INFILE1" --right "$INFILE2" -c "$CPU" -o "$WORKDIR/$BASE"
        else
            echo "Trimmed reads already exist; skipping."
        fi
        check_file "$LEFTTRIM" "Trimmed left reads"
        check_file "$RIGHTTRIM" "Trimmed right reads"

        echo
        echo "--- Step 2: Filtering reads ---"
        if [[ ! -s "$LEFT" || ! -s "$RIGHT" ]]; then
            AAFTF filter -c "$CPU" -o "$WORKDIR/$BASE" --left "$LEFTTRIM" --right "$RIGHTTRIM" --aligner bbduk
        else
            echo "Filtered reads already exist; skipping."
        fi
        check_file "$LEFT" "Filtered left reads"
        check_file "$RIGHT" "Filtered right reads"

        echo
        echo "--- Finished preprocessing for $DES ---"
    } > "$log_temp" 2>&1

    sed -E "s/\x1B\[[0-9;]*[mK]//g" "$log_temp" > "$log_file"
    rm -f "$log_temp"
done < <(sed 's/\r$//' "$INPUT_FILE")

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
    echo "  --maxMemory <GB>        Max memory in GB (auto-detected if not provided)"
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
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' does not exist."
    exit 1
fi

# --- Check dependencies ---
for cmd in funannotate; do
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
RNA_DIR=$(sed -n 's/^#Location of RNAseq Reads=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' "$INPUT_FILE" | head -n1 | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')

if [[ -z "$RNA_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Could not extract required directories from $INPUT_FILE"
    exit 1
fi

# --- Validate input directory ---
[[ "$RNA_DIR" != "NA" ]] && [[ ! -d "$RNA_DIR" ]] && echo "Error: $RNA_DIR does not exist." && exit 1

# --- Helper function ---
check_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -s "$file" ]]; then
        echo "Error: $desc ($file) not found or empty."
        exit 1
    fi
}

# --- Directory setup ---
LOGDIR="$OUTPUT_DIR/logs"
ASM_DIR="$OUTPUT_DIR/assembled_genomes"
ANNO_DIR="$OUTPUT_DIR/funannotate_genomes"

mkdir -vp "$ANNO_DIR" "$LOGDIR" 

# --- Preflight check: sanitized DES collisions ---
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

    RNAseq_FLAG=0
    [[ "$RNA1" != "NA" ]] && RNAseq_FLAG=1

    # --- Define file paths #
    BASE=${DES//[^a-zA-Z0-9._-]/_}

    POLISHDIR="$ASM_DIR/$BASE.$ASSEMBLER.polishout"
    SORTED="$POLISHDIR/$BASE.sorted.fasta"

    ANNO_SAMPLE_DIR="$ANNO_DIR/$BASE.$ASSEMBLER"
    mkdir -vp "$ANNO_SAMPLE_DIR" 
    cd "$ANNO_SAMPLE_DIR"

    CLEAN="$ANNO_SAMPLE_DIR/$BASE.funannotate01.fasta"
    SORT2="$ANNO_SAMPLE_DIR/$BASE.funannotate02.fasta"
    MASK="$ANNO_SAMPLE_DIR/$BASE.funannotate03.fasta"
    PREDICT="$ANNO_SAMPLE_DIR/$BASE.funannotate04"
    PREDICT_STATS="$PREDICT/predict_results/${SPECIES}_$BASE.stats.json"
    ANNO_OUT="$PREDICT/annotate_results/${SPECIES}_${BASE}.stats.json"

    timestamp=$(date +%Y%m%d_%H%M%S)
    log_temp="$LOGDIR/part4.$BASE.$timestamp.log.tmp"
    log_file="$LOGDIR/part4.$BASE.$timestamp.log"

    {
        echo "--- Annotation for $DES ($ASSEMBLER) ---"
        echo "Input file: $INPUT_FILE | Threads: $THREADS | Memory: $MAX_MEMORY GB"
        echo "Started at: $(date)"
        check_file "$SORTED" "Sorted assembly input"

        echo
        echo "--- Step 1: funannotate clean ---"
        if [[ ! -s "$CLEAN" ]]; then
            funannotate clean -i "$SORTED" -o "$CLEAN"
        else
            echo "$CLEAN exists; skipping clean step."
        fi
        check_file "$CLEAN" "Cleaned scaffolds"

        echo
        echo "--- Step 2: funannotate sort ---"
        if [[ ! -s "$SORT2" ]]; then
            funannotate sort -i "$CLEAN" -o "$SORT2" --minlen 50 --base scaffold
        else
            echo "$SORT2 exists; skipping sort step."
        fi
        check_file "$SORT2" "Sorted scaffolds"

        echo
        echo "--- Step 3: funannotate mask ---"
        if [[ ! -s "$MASK" ]]; then
            funannotate mask -i "$SORT2" -o "$MASK" --cpus "$THREADS"
        else
            echo "$MASK exists; skipping mask step."
        fi
        check_file "$MASK" "Masked scaffolds"
        
        if [[ "$RNAseq_FLAG" -eq 1 ]]; then        
            echo
            echo "--- Step 4: funannotate train ---"
            if [[ ! -s "$PREDICT/training/funannotate_train.pasa.gff3" ]]; then
                if [[ "$RNA2" != "NA" ]]; then
                    funannotate train -i "$MASK" \
                        -o "$PREDICT" \
                        --left "${RNA_DIR}/${RNA1}" \
                        --right "${RNA_DIR}/${RNA2}" \
                        --species "$SPECIES" \
                        --strain "$BASE" \
                        --cpus "$THREADS" \
                        --pasa_db sqlite
                else
                    funannotate train -i "$MASK" \
                        -o "$PREDICT" \
                        --single "${RNA_DIR}/${RNA1}" \
                        --species "$SPECIES" \
                        --strain "$BASE" \
                        --cpus "$THREADS"
                fi
            else
                echo "Training output seems to be present and complete; (e.g the file $PREDICT/training/funannotate_train.pasa.gff3 has been produced); skipping training step."
            fi 
        else
            echo "No RNA-seq data detected for $DES; skipping training step."
        fi

        echo
        echo "--- Step 5: funannotate predict ---"
        if [[ ! -f "$PREDICT_STATS" ]]; then
            funannotate predict -i "$MASK" -o "$PREDICT" \
                --species "$SPECIES" --strain "$BASE" \
                --cpus "$THREADS" --name "$BASE"
        else
            echo "Prediction output exists; skipping predict step."
        fi
        check_file "$PREDICT_STATS" "Prediction stats JSON"

        echo
        echo "--- Step 5.5: eggnog mapper ---"
        if [[ ! -f "$PREDICT/annotate_misc/eggnog.emapper.annotations" ]]; then
            eggnog_only.py -i "$PREDICT/predict_results/${SPECIES}_${BASE}.proteins.fa" -o "$PREDICT/annotate_misc/eggnog" --cpus "$THREADS"
        else
            echo "Eggnog mapper output exists; skipping eggnog mapper."
        fi
        check_file "$PREDICT/annotate_misc/eggnog.emapper.annotations" "Eggnog mapper annotations output"

        echo "--- Step 6: funannotate annotate ---"
        if [[ ! -f "$ANNO_OUT" ]]; then
            funannotate annotate -i "$PREDICT" -o "$PREDICT" \
                --species "$SPECIES" --strain "$BASE" \
                --cpus "$THREADS"
        else
            echo "Annotation output exists; skipping annotate step."
        fi
        check_file "$ANNO_OUT" "Annotation stats JSON"

        echo
        echo "--- Completed annotation for $DES ---"
    } > "$log_temp" 2>&1

    sed -E "s/\x1B\[[0-9;]*[mK]//g" "$log_temp" > "$log_file"
    rm -f "$log_temp"

done < "$INPUT_FILE"
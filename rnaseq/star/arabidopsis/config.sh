#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
# set -eou pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="arabidopsis_PRJNA1225620"
export SPECIES="Arabidopsis thaliana"

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ========== 2. BASIC PROJECT PATHS ==========
export ROOT_DIR="${SCRIPT_DIR}"
export DATA_DIR="${ROOT_DIR}/data"
export REF_DIR="${ROOT_DIR}/reference"
export SCRIPTS_DIR="${ROOT_DIR}/scripts"
export TRIM_DIR="${ROOT_DIR}/data_trimmed"

export BASE_RESULTS_DIR="${ROOT_DIR}/results"
export OUT_DIR="${BASE_RESULTS_DIR}/${PROJECT_NAME}"

# STAR specific directories
export STAR_INDEX_DIR="${REF_DIR}/star_index"
export ALIGN_DIR="${OUT_DIR}/star_alignments"
export COUNTS_DIR="${OUT_DIR}/counts"
export QC_DIR="${OUT_DIR}/fastqc"

mkdir -p "${DATA_DIR}" "${REF_DIR}" "${STAR_INDEX_DIR}" "${SCRIPTS_DIR}" "${TRIM_DIR}" "${ALIGN_DIR}" "${COUNTS_DIR}" "${QC_DIR}" "${BASE_RESULTS_DIR}" "${OUT_DIR}"

# ========== 3. DATA ACQUISITION URLs & LISTS ==========
export SRR_LIST="${ROOT_DIR}/SRR_acc_list.txt"

# STAR requires the DNA genome and the GTF Annotation
export GENOME_URL="http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
export GENOME_FA="${REF_DIR}/athal_genome.fa"

export ANNOTATION_URL="http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz"
export ANNOTATION_GTF="${REF_DIR}/athal_annotation.gtf"

# ========== 4. SAMPLE METADATA ==========
export METADATA_FILE="${ROOT_DIR}/metadata.csv"

SAMPLES=()
TREATMENTS=()
while IFS=',' read -r samp treat || [ -n "$samp" ]; do
    samp=$(echo "$samp" | tr -d '\r' | tr -d ' ')
    treat=$(echo "$treat" | tr -d '\r' | tr -d ' ')
    if [[ -n "$samp" && "$samp" != "Sample" ]]; then
        SAMPLES+=("$samp")
        TREATMENTS+=("$treat")
    fi
done < "${METADATA_FILE}"

export SAMPLE_LIST=$(IFS=,; echo "${SAMPLES[*]}")
export TREATMENT_LIST=$(IFS=,; echo "${TREATMENTS[*]}")
export CONTROL_TREATMENT="WT"

# ========== 5. QUANTIFICATION PARAMETERS ==========
export THREADS=8
export READ_LENGTH=150  # Change this if your read length is 100 or 50
export OVERHANG=$((READ_LENGTH - 1))

echo "------------------------------------------------"
echo "  STAR RNA-seq Pipeline Config Loaded           "
echo "  Project/Run Name: ${PROJECT_NAME}             "
echo "  Output Directory: ${OUT_DIR}                  "
echo "------------------------------------------------"

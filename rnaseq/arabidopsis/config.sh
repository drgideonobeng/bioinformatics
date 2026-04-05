#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
# set -eo pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="Arabidopsis_Salmon_RNAseq"
export SPECIES="Arabidopsis thaliana"

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ========== 2. BASIC PROJECT PATHS ==========
export ROOT_DIR="${SCRIPT_DIR}"

# Input/Output Directories
export DATA_DIR="${ROOT_DIR}/data"
export REF_DIR="${ROOT_DIR}/reference"
export INDEX_DIR="${REF_DIR}/salmon_index"
export RESULTS_DIR="${ROOT_DIR}/results/quants"
export SCRIPTS_DIR="${ROOT_DIR}/scripts"

# Main OUT_DIR so R knows where to save the .rds files and plots
export OUT_DIR="${ROOT_DIR}/results"
export PLOT_DIR="${OUT_DIR}/plots"

# Ensure all result and reference directories exist safely
mkdir -p "${DATA_DIR}" "${REF_DIR}" "${INDEX_DIR}" "${RESULTS_DIR}" "${SCRIPTS_DIR}" "${PLOT_DIR}"

# ========== 3. DATA ACQUISITION URLs & LISTS ==========
# Pointing directly to your SRA Accession List
export SRR_LIST="${ROOT_DIR}/SRR_acc_list.txt"

# Updated to the current Ensembl TAIR10.1 reference files
export TRANSCRIPTOME_URL="http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
export TRANSCRIPTOME_FA="${REF_DIR}/athal.fa.gz"

export GENOME_URL="http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
export GENOME_FA="${REF_DIR}/athal_genome.fa.gz"

# ========== 4. SAMPLE METADATA (Cross-Language) ==========
export METADATA_FILE="${ROOT_DIR}/metadata.csv"

# Initialize empty arrays
SAMPLES=()
TREATMENTS=()

# Read the CSV line by line (Only needed for your R metadata pipeline now)
while IFS=',' read -r samp treat || [ -n "$samp" ]; do
    # 1. Scrub hidden Windows carriage returns (\r) and accidental spaces
    samp=$(echo "$samp" | tr -d '\r' | tr -d ' ')
    treat=$(echo "$treat" | tr -d '\r' | tr -d ' ')

    # 2. Skip the header row and any completely blank lines
    if [[ -n "$samp" && "$samp" != "Sample" ]]; then
        SAMPLES+=("$samp")
        TREATMENTS+=("$treat")
    fi
done < "${METADATA_FILE}"

# Export the comma-separated strings so your R scripts can still read them perfectly
export SAMPLE_LIST=$(IFS=,; echo "${SAMPLES[*]}")
export TREATMENT_LIST=$(IFS=,; echo "${TREATMENTS[*]}")

export CONTROL_TREATMENT="Water"

# ========== 5. QUANTIFICATION PARAMETERS ==========
export THREADS=8
export LIB_TYPE="A"  # 'A' tells Salmon to automatically infer library type
export KMER_SIZE=31  # Default k-mer size for Salmon index

echo "------------------------------------------------"
echo "  Salmon RNA-seq Pipeline Config Loaded         "
echo "  Project: ${PROJECT_NAME}                      "
echo "  Species: ${SPECIES}                           "
echo "  Root: ${ROOT_DIR}                             "
echo "------------------------------------------------"

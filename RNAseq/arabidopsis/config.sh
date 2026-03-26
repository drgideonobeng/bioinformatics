#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
set -eo pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="Arabidopsis_Salmon_RNAseq"
export SPECIES="Arabidopsis thaliana"

# ========== 2. BASIC PROJECT PATHS ==========
# ROOT_DIR is the directory where this config.sh file lives
export ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input/Output Directories
export DATA_DIR="${ROOT_DIR}/data"
export REF_DIR="${ROOT_DIR}/reference"
export INDEX_DIR="${REF_DIR}/salmon_index"
export RESULTS_DIR="${ROOT_DIR}/results/quants"
export SCRIPTS_DIR="${ROOT_DIR}/scripts"

# Ensure all result and reference directories exist safely
mkdir -p "${DATA_DIR}" "${REF_DIR}" "${INDEX_DIR}" "${RESULTS_DIR}" "${SCRIPTS_DIR}"

# ========== 3. DATA ACQUISITION URLs & SAMPLES ==========
# Reference Transcriptome (Ensembl)
export TRANSCRIPTOME_URL="ftp://ftp.ensemblgenomes.org/pub/plants/release-28/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.28.cdna.all.fa.gz"
export TRANSCRIPTOME_FA="${REF_DIR}/athal.fa.gz"

# Sample IDs (Subset of PRJDB2508 for demonstration)
export SAMPLES=("DRR016125" "DRR016126" "DRR016127" "DRR016128")
export SRA_BASE_URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR016"

# ========== 4. QUANTIFICATION PARAMETERS ==========
export THREADS=4
export LIB_TYPE="A"  # 'A' tells Salmon to automatically infer library type (stranded/unstranded)
export KMER_SIZE=31  # Default k-mer size for Salmon index

echo "------------------------------------------------"
echo "  Salmon RNA-seq Pipeline Config Loaded         "
echo "  Project: ${PROJECT_NAME}                      "
echo "  Species: ${SPECIES}                           "
echo "  Samples: ${#SAMPLES[@]} queued for processing "
echo "  Root: ${ROOT_DIR}                             "
echo "------------------------------------------------"

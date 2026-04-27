#!/usr/bin/env bash
set -eou pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="melanoma_crispr_screen"
# Get absolute path to the directory this script lives in
export ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ========== 2. PROJECT DIRECTORIES ==========
export DATA_DIR="${ROOT_DIR}/data"
export RAW_DIR="${DATA_DIR}/raw"
export TRIM_DIR="${DATA_DIR}/trimmed"
export REF_DIR="${DATA_DIR}/reference"
export SCRIPTS_DIR="${ROOT_DIR}/scripts"

export OUT_DIR="${ROOT_DIR}/results/${PROJECT_NAME}"
export QC_DIR="${OUT_DIR}/01_qc"
export COUNT_DIR="${OUT_DIR}/02_counts"
export TEST_DIR="${OUT_DIR}/03_testing"

# Ensure directories exist
mkdir -p "${RAW_DIR}" "${TRIM_DIR}" "${REF_DIR}" "${QC_DIR}" "${COUNT_DIR}" "${TEST_DIR}"

# ========== 3. DATASET URLs (Fujihara et al. subset for demonstration) ==========
export T0_CTRL_URL="https://zenodo.org/records/5750854/files/T0-Control.fastq.gz"
export T8_VEH_URL="https://zenodo.org/records/5750854/files/T8-Vehicle.fastq.gz"
export T8_DRUG_URL="https://zenodo.org/records/5750854/files/T8-APR-246.fastq.gz"
export LIBRARY_URL="https://zenodo.org/records/5750854/files/Brunello_library.csv"

# ========== 4. EXPERIMENT METADATA ==========
# Adapter sequence to trim (5' adapter for the lentiGuide-Puro vector)
export ADAPTER_5P="TCTTGTGGAAAGGACGAAACACCG" 

# MAGeCK labels
export CTRL_SAMPLES="T8-Vehicle"
export TREAT_SAMPLES="T8-APR-246"
export LIBRARY_FILE="${REF_DIR}/Brunello_library.csv"

export THREADS=8

echo "------------------------------------------------"
echo "  CRISPR Pipeline Config Loaded                 "
echo "  Project Name: ${PROJECT_NAME}                 "
echo "  Output Directory: ${OUT_DIR}                  "
echo "------------------------------------------------"

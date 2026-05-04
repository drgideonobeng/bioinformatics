#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
# set -eo pipefail # Activated your safety net!

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="Mesp1_E80"
export GEO_ACCESSION="GSM4816083"
export GEO_PREFIX="${GEO_ACCESSION}_${PROJECT_NAME}"

# ========== 2. BASIC PROJECT PATHS ==========
# ROOT_DIR is the directory where this config.sh file lives
export ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input/Output Directories
export DATA_MATRIX_DIR="${ROOT_DIR}/data/filtered_feature_bc_matrix"
export OBJ_DIR="${ROOT_DIR}/results/objects"
export PLOT_DIR="${ROOT_DIR}/results/plots"
export TABLE_DIR="${ROOT_DIR}/results/tables"

# Ensure all result directories exist safely via Bash
mkdir -p "${DATA_MATRIX_DIR}" "${OBJ_DIR}" "${PLOT_DIR}" "${TABLE_DIR}"

# ========== 3. QC & ANALYSIS THRESHOLDS ==========
export MIN_CELLS=3
export MIN_FEATURES_BASE=200
export MT_PATTERN="^mt-"

export MIN_GENES=2000
export MAX_GENES=6000
export MAX_MT_PERCENT=3

export PCA_DIMS=20
export CLUSTER_RES=0.4
export UMAP_TITLE="${PROJECT_NAME} lineage cells"

# ========== 4. COMPUTATIONAL RESOURCES ==========
export THREADS=8

echo "------------------------------------------------"
echo "  Universal scRNA-seq Pipeline Config Loaded    "
echo "  Project: ${PROJECT_NAME}                      "
echo "  Root: ${ROOT_DIR}                             "
echo "  Resolution: ${CLUSTER_RES}                    "
echo "------------------------------------------------"

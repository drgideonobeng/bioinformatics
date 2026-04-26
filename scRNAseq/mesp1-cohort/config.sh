#!/usr/bin/env bash
# set -eo pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="Mesp1_Cohort" 

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

# ========== 3. GLOBAL QC & ANALYSIS THRESHOLDS ==========
export MIN_CELLS=3
export MIN_FEATURES_BASE=200
export MT_PATTERN="^mt-"

export MIN_GENES=200
export MAX_GENES=5500
export MAX_MT_PERCENT=5

export PCA_DIMS=20
export CLUSTER_RES=0.5
export UMAP_TITLE="${PROJECT_NAME} Integrated Lineage"

export CLUSTER_NAME="Endothelial"

# ========== 4. COMPUTATIONAL RESOURCES ==========
export THREADS=8

echo "------------------------------------------------"
echo "  Multi-Sample Pipeline Config Loaded           "
echo "  Project: ${PROJECT_NAME}                      "
echo "  Root: ${ROOT_DIR}                             "
echo "------------------------------------------------"

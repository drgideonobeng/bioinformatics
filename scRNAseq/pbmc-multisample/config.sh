#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
# set -eo pipefail

# ========== 1. PROJECT IDENTIFICATION ==========
export PROJECT_NAME="PBMC_Multi_Sample"
export MT_PATTERN="^MT-" # Use ^MT- for human, ^mt- for mouse

# ========== 2. BASIC PROJECT PATHS ==========
# ROOT_DIR is the directory where this config.sh file lives
export ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input/Output Directories 
export DATA_DIR_UNTREATED="${ROOT_DIR}/data/untreated"
export DATA_DIR_TREATED="${ROOT_DIR}/data/treated"

export OBJ_DIR="${ROOT_DIR}/results/objects"
export PLOT_DIR="${ROOT_DIR}/results/plots"
export TABLE_DIR="${ROOT_DIR}/results/tables"

# Ensure all result directories exist safely via Bash
mkdir -p "${OBJ_DIR}" "${PLOT_DIR}" "${TABLE_DIR}"

# ========== 3. QC & ANALYSIS THRESHOLDS ==========
# Seurat Object Initialization Prefilters
export MIN_CELLS=3
export MIN_FEATURES_BASE=200

# Rigorous Cell Filtering
export MIN_GENES=200
export MAX_GENES=5000  
export MAX_MT_PERCENT=15

# Clustering and DimRed parameters
export PCA_DIMS=20
export CLUSTER_RES=0.5
export UMAP_TITLE="Untreated vs Treated PBMCs"
export TOP_GENES=2

# ========== 4. COMPUTATIONAL RESOURCES ==========
export THREADS=4

# ========== 5. DIFFERENTIAL EXPRESSION METADATA ==========
# These are the exact labels Seurat will use for A/B testing in Step 09.
# NOTE: Make sure these match the EXACT case of your folder names!
export CONTROL_GROUP_NAME="untreated"
export TEST_GROUP_NAME="treated"

echo "------------------------------------------------"
echo "  Universal scRNA-seq Pipeline Config Loaded    "
echo "  Project: ${PROJECT_NAME}                      "
echo "  Comparing: ${CONTROL_GROUP_NAME} vs ${TEST_GROUP_NAME} "
echo "  Root: ${ROOT_DIR}                             "
echo "------------------------------------------------"

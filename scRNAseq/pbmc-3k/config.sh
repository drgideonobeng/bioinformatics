#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
set -euo pipefail

# ========== BASIC PROJECT PATHS ==========
# ROOT_DIR is the directory where this config.sh file lives
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input/Output Directories (Updated to match 10x tarball extraction)
export DATA_DIR="${ROOT_DIR}/data/filtered_gene_bc_matrices/hg19"
export OBJ_DIR="${ROOT_DIR}/results/objects"
export PLOT_DIR="${ROOT_DIR}/results/plots"

# Reference Directory
export REF_DIR="${ROOT_DIR}/ref"

# Ensure all result directories exist
mkdir -p "${OBJ_DIR}" "${PLOT_DIR}" "${REF_DIR}"

# ========== QC & ANALYSIS THRESHOLDS ==========
# These will be pulled into R using Sys.getenv()
export MIN_GENES=200
export MAX_GENES=2500
export MAX_MT_PERCENT=5

# Clustering and DimRed parameters
export PCA_DIMS=10
export CLUSTER_RES=0.5

# ========== COMPUTATIONAL RESOURCES ==========
export THREADS=4

echo "------------------------------------------------"
echo "  PBMC scRNA-seq Pipeline Config Loaded         "
echo "  Root: ${ROOT_DIR}                             "
echo "  Data: ${DATA_DIR}                             "
echo "------------------------------------------------"

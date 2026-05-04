#!/bin/bash
# ==============================================================================
# config.sh
# Master configuration file for the scATAC-seq modular pipeline
# ==============================================================================

# 1. Project & Data Source
export PROJECT_NAME="PBMC_10k"

# ========== 2. BASIC PROJECT PATHS ==========
# ROOT_DIR is the directory where this config.sh file lives
export ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Input/Output Directories
export DATA_MATRIX_DIR="${ROOT_DIR}/data/raw_matrix"
export OBJ_DIR="${ROOT_DIR}/results/objects"
export PLOT_DIR="${ROOT_DIR}/results/plots"
export TABLE_DIR="${ROOT_DIR}/results/tables"

# 3. Filtering Parameters (scATAC-seq specific)
export MIN_PEAKS="3000"
export MAX_PEAKS="100000"
export MIN_TSS_ENRICHMENT="2"
export MAX_NUCLEOSOME_SIGNAL="4"

# 4. Dimensionality Reduction & Clustering
export LSI_DIMS_START="2" # Skip 1 in ATAC-seq (correlates with depth)
export LSI_DIMS_END="30"
export CLUSTER_RES="0.8"

# 5. ATAC Object Parameters
export MIN_CELLS="10"
export MIN_FEATURES="200"
export GENOME_BUILD="hg38"

# 6. Data Download URLs
export URL_MATRIX="https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
export URL_FRAGMENTS="https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz"
export URL_FRAG_INDEX="https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"

echo "=================================================="
echo "Configuration loaded for project: $PROJECT_NAME"
echo "Root Directory : $ROOT_DIR"
echo "Minimum TSS Enrichment : $MIN_TSS_ENRICHMENT"
echo "Maximum Nucleosome Signal : $MAX_NUCLEOSOME_SIGNAL"
echo "=================================================="

#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
set -euo pipefail

echo "=================================================="
echo "      Starting PBMC scRNA-seq Pipeline            "
echo "=================================================="

# 1. Load the environment variables
# This guarantees that all R scripts inherit these settings
if [ -f "config.sh" ]; then
    source config.sh
else
    echo "ERROR: config.sh not found. Please run this script from the project root."
    exit 1
fi

# 2. Execute the R scripts in order
echo "[1/9] Checking for raw data..."
Rscript scripts/01_data_download.R

echo "[2/9] Creating Seurat object..."
Rscript scripts/02_create_seurat_obj.R

echo "[3/9] Generating pre-filter QC plots..."
Rscript scripts/03_qc_visualize.R

echo "[4/9] Filtering and normalizing cells..."
Rscript scripts/04_filter_cells.R

echo "[5/9] Normalizing and scaling data..."
Rscript scripts/05_normalize_data.R

echo "[6/9] Running PCA and UMAP clustering..."
Rscript scripts/05_run_dim_reduction.R

echo "[7/9] Generating clustering plots..."
Rscript scripts/06_visualize_clusters.R

echo "[8/9] Finding cluster marker genes..."
Rscript scripts/07_find_markers.R

echo "[9/9] Annotating final cell types..."
Rscript scripts/08_annotate_cells.R

echo "=================================================="
echo "  PIPELINE COMPLETE! Check the results/ folder.   "
echo "=================================================="

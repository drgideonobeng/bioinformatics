#!/bin/bash
# ==============================================================================
# run_pipeline.sh
# Executes the scATAC-seq R scripts in sequence
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e 

# 1. Source the configuration file to load environment variables
source config.sh

echo "========================================="
echo " Starting scATAC-seq Pipeline"
echo "========================================="

# 2. Execute R scripts in sequence
# (Assuming your R scripts are in a folder named 'scripts' and are executable)

echo "Step 1: Downloading data..."
Rscript scripts/01_atac_data_download.R

echo "Step 2: Creating Seurat/Signac Object..."
Rscript scripts/02_create_atac_obj.R

echo "Step 3: Calculating QC Metrics..."
Rscript scripts/03_atac_qc_visualize.R

echo "Step 4: Filtering Cells..."
Rscript scripts/04_filter_atac_cells.R

echo "Step 5: Normalization (TF-IDF) & SVD..."
Rscript scripts/05_normalize_lsi.R

echo "Step 6: UMAP & Clustering..."
Rscript scripts/06_atac_umap_cluster.R

echo "Step 7: Calculating Gene Activity..."
Rscript scripts/07_gene_activity.R

echo "========================================="
echo " Pipeline Complete! Check the '$PLOT_DIR' folder."
echo "========================================="

#!/usr/bin/env bash

# Exit immediately if any command fails
set -e

echo "======================================================"
echo "  Starting Single-Cell RNA-Seq Analysis Pipeline...   "
echo "======================================================"

# 1. Load the configuration variables
source ./config.sh

echo ""
echo ">>> Part 0: Data Fetching <<<"
Rscript scripts/00_fetch_benchmark_data.R

echo ""
echo ">>> Part 1: Pre-Processing & Quality Control <<<"
Rscript scripts/01_load_and_merge_data.R
Rscript scripts/02_quality_control.R
Rscript scripts/03_normalize.R

echo ""
echo ">>> Part 2: Dimensionality Reduction & Clustering <<<"
Rscript scripts/04_scale_and_pca.R
Rscript scripts/05_cluster_and_umap.R
Rscript scripts/06_find_markers.R
Rscript scripts/07_visualize_markers.R

echo ""
echo ">>> Part 3: Annotation & Differential Expression <<<"
# Running both manual and automated annotations for completeness
Rscript scripts/08a_annotate_clusters_manual.R
Rscript scripts/08b_annotate_clusters_automated.R
Rscript scripts/09_differential_expression.R
Rscript scripts/10_volcano_plot.R
Rscript scripts/11_pathway_analysis.R

echo ""
echo "======================================================"
echo " PIPELINE COMPLETE"
echo " All objects saved to:   ${OBJ_DIR}"
echo " All plots saved to:     ${PLOT_DIR}"
echo " All tables saved to:    ${TABLE_DIR}"
echo "======================================================"

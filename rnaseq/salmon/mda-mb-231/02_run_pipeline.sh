#!/usr/bin/env bash

#!/usr/bin/env bash
# set -e: exit on error; -u: exit on unset variables; -o pipefail: catch errors in pipes
set -euo pipefail

# 1. SOURCE THE CONFIG! This passes all your dynamic paths to R
source config.sh

echo "########################################################"
echo " INITIALIZING AUTOMATED R DOWNSTREAM PIPELINE"
echo " Project: ${PROJECT_NAME}"
echo " Output Dir: ${OUT_DIR}"
echo "########################################################"

echo "=> [1/7] Importing Data & Setup (tximport)..."
Rscript rscripts/01_import_setup.R

echo "=> [2/7] Executing DESeq2 Model..."
Rscript rscripts/02_model_execution.R

echo "=> [3/7] Generating QC Visualizations..."
Rscript rscripts/03_qc_visualizations.R

echo "=> [4/7] Extracting Differential Expression Results..."
Rscript rscripts/04_extract_results.R

echo "=> [5/7] Annotating Genes..."
# Notice: I removed the hardcoded CLE9 paths here! We will make script 05 smart enough to find the tables dynamically.
Rscript rscripts/05_annotate_genes.R

echo "=> [6/7] Running GO Enrichment..."
Rscript rscripts/06_go_enrichment.R

echo "=> [7/7] Running Gene Set Enrichment Analysis (GSEA)..."
Rscript rscripts/07_run_gsea.R

echo "########################################################"
echo " R PIPELINE COMPLETE!"
echo " All statistics, annotations, and plots are in: ${OUT_DIR}"
echo "########################################################"

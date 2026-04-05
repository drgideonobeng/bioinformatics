#!/usr/bin/env bash

# Exit immediately if any command fails
set -e

echo "########################################################"
echo " INITIALIZING AUTOMATED R DOWNSTREAM PIPELINE"
echo "########################################################"

# Ensure the results directory variable is set so Step 1 knows where to look
export RESULTS_DIR="results/quants"

echo "=> [1/7] Importing Data & Setup (tximport)..."
Rscript rscripts/01_import_setup.R

echo "=> [2/7] Executing DESeq2 Model..."
Rscript rscripts/02_model_execution.R

echo "=> [3/7] Generating QC Visualizations..."
Rscript rscripts/03_qc_visualizations.R

echo "=> [4/7] Extracting Differential Expression Results..."
Rscript rscripts/04_extract_results.R

echo "=> [5/7] Annotating Genes..."
# Run the annotation script on both the full genome list and the significant DEGs
Rscript rscripts/05_annotate_genes.R results/DE_results/CLE9_vs_Water_all_genes.csv
Rscript rscripts/05_annotate_genes.R results/DE_results/CLE9_vs_Water_sig_DEGs.csv

echo "=> [6/7] Running GO Enrichment..."
Rscript rscripts/06_go_enrichment.R

echo "=> [7/7] Running Gene Set Enrichment Analysis (GSEA)..."
Rscript rscripts/07_run_gsea.R

echo "########################################################"
echo " R PIPELINE COMPLETE!"
echo " All statistics, annotations, and plots are in the results/ folder."
echo "########################################################"

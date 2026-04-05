#!/bin/bash
# ==========================================
# Project: 03-yeast-metabolism (Glucose vs Galactose)
# ==========================================

PROJECT_DIR=$(pwd)
DATA_DIR="$PROJECT_DIR/data"
REF_DIR="$DATA_DIR/reference"
RESULTS_DIR="$PROJECT_DIR/results"

# Metabolic switch samples
SAMPLES=("SRR1543118" "SRR1543119")

THREADS=4

# Reference Files
TRANSCRIPTOME_FA="$REF_DIR/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa"
SALMON_INDEX="$REF_DIR/salmon_index"
GENOME_GTF="$REF_DIR/Saccharomyces_cerevisiae.R64-1-1.104.gtf"

mkdir -p "$DATA_DIR/raw" "$REF_DIR" "$RESULTS_DIR/salmon_quants"

echo "----------------------------------------------------------------"
echo "Config Loaded: Yeast Metabolic Switch Project"
echo "Comparing: Glucose (Control) vs Galactose (Treated)"
echo "----------------------------------------------------------------"

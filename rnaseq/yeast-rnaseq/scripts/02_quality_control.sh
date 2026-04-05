#!/bin/bash
source config.sh

# Define directories
RAW_DIR="$DATA_DIR/raw"
QC_OUT="$RESULTS_DIR/fastqc"
MQC_OUT="$RESULTS_DIR/multiqc"

# Create output directories
mkdir -p "$QC_OUT" "$MQC_OUT"

echo "=== Running FastQC on Raw Reads ==="
# -t uses multiple CPU cores for speed
fastqc -t $THREADS "$RAW_DIR"/*.fastq.gz -o "$QC_OUT"

echo "=== Aggregating Reports with MultiQC ==="
# multiqc scans the directory for any QC logs and builds a report
multiqc "$QC_OUT" -o "$MQC_OUT"

echo "----------------------------------------------------------------"
echo "[OK] QC Complete."
echo "View your report here: $MQC_OUT/multiqc_report.html"
echo "----------------------------------------------------------------"

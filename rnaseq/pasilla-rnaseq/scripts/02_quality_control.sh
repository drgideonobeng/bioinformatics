#!/bin/bash

RAW_DIR="data/raw"
OUT_DIR="results/fastqc"

echo "=== Running Quality Control (FastQC) ==="

# Run FastQC on all compressed fastq files
for file in ${RAW_DIR}/*.fastq.gz; do
    echo "Processing: $file"
    fastqc "$file" -o "$OUT_DIR"
done

echo "[OK] HTML reports generated in ${OUT_DIR}"

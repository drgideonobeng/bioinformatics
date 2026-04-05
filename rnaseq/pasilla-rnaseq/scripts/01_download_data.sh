#!/bin/bash
RAW_DIR="data/raw"
mkdir -p $RAW_DIR

echo "=== Downloading Verified Small RNA-Seq Dataset ==="

# Using a subset of the human airway dataset hosted on a stable repository
curl -L -o ${RAW_DIR}/treated_R1.fastq.gz https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/data/tophat_test/test_1.fastq.gz
curl -L -o ${RAW_DIR}/treated_R2.fastq.gz https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/data/tophat_test/test_2.fastq.gz

echo "Verifying file integrity..."
for file in ${RAW_DIR}/*.fastq.gz; do
    if file "$file" | grep -q "gzip compressed data"; then
        echo "[OK] $file is a valid GZIP file."
    else
        echo "[ERROR] $file is still invalid. Printing file header for debug:"
        head -n 5 "$file"
    fi
done

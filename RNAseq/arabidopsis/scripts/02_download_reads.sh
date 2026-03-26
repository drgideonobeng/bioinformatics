#!/usr/bin/env bash

# Load the master configuration file
source "$(dirname "$0")/../config.sh"

echo "======================================================"
echo " Step 2: Downloading RNA-seq Reads"
echo "======================================================"

for SAMP in "${SAMPLES[@]}"; do
    echo "=> Fetching reads for sample: ${SAMP}"
    
    # Create a specific directory for this sample's raw data
    mkdir -p "${DATA_DIR}/${SAMP}"
    
    # Download Forward (_1) and Reverse (_2) reads
    curl -L -f -C - "${SRA_BASE_URL}/${SAMP}/${SAMP}_1.fastq.gz" -o "${DATA_DIR}/${SAMP}/${SAMP}_1.fastq.gz"
    curl -L -f -C - "${SRA_BASE_URL}/${SAMP}/${SAMP}_2.fastq.gz" -o "${DATA_DIR}/${SAMP}/${SAMP}_2.fastq.gz"
done

echo "=> All raw data downloaded successfully to: ${DATA_DIR}"

#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 2: Downloading RNA-seq Reads via SRA Toolkit"
echo "======================================================"

# Check if the SRR list exists
if [ ! -f "${SRR_LIST}" ]; then
    echo "ERROR: Could not find ${SRR_LIST}!"
    exit 1
fi

echo "=> Processing SRA Accessions..."

# Loop through each SRR number in your text file
for SRR in $(cat "${SRR_LIST}"); do
    echo "------------------------------------------------"
    echo " Processing sample: ${SRR}"
    
    # Create a clean subfolder for each sample
    SAMP_DIR="${DATA_DIR}/${SRR}"
    mkdir -p "${SAMP_DIR}"
    
    # Check if the final compressed files already exist
    # If they do, 'continue' skips to the NEXT sample in the list!
    if [ -f "${SAMP_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${SAMP_DIR}/${SRR}_2.fastq.gz" ]; then
        echo "   -> Final .fastq.gz files already exist. Skipping ${SRR} entirely."
        continue
    fi

    # --- 1. PREFETCH DATA ---
    echo "   -> Safely prefetching compressed SRA files..."
    # Now this only runs if the final fastq files don't exist
    prefetch "${SRR}" --output-directory "${DATA_DIR}" 

    # --- 2. EXTRACT & COMPRESS DATA ---
    echo "   -> Extracting to FASTQ format..."
    fasterq-dump "${SRR}" --split-files --threads ${THREADS} --outdir "${SAMP_DIR}"
    
    echo "   -> Compressing files to .fastq.gz to save space for Salmon..."
    gzip "${SAMP_DIR}/${SRR}_1.fastq"
    gzip "${SAMP_DIR}/${SRR}_2.fastq"
    
    echo "   -> ${SRR} complete!"
done

echo "======================================================"
echo "=> All SRA data processed successfully in: ${DATA_DIR}"
echo "======================================================"

#!/bin/bash
source config.sh

echo "=== Downloading Metabolic Switch Samples ==="

for ACC in "${SAMPLES[@]}"; do
    echo "Processing $ACC..."
    # Download from NCBI
    prefetch "$ACC"
    # Convert to FastQ (split into R1 and R2)
    fasterq-dump --split-files "$ACC" -O "$DATA_DIR/raw"
    # Compress to save space (Salmon reads .gz)
    pigz "$DATA_DIR/raw/${ACC}"_*.fastq
done

echo "[OK] Raw data downloaded and compressed."

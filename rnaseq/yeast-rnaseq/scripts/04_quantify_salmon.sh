#!/bin/bash
source config.sh

RAW_DIR="$DATA_DIR/raw"
OUT_DIR="$RESULTS_DIR/salmon_quants"

echo "=== Running Salmon Quantification (Single-End) ==="

for ACC in "${SAMPLES[@]}"; do
    echo "Quantifying $ACC..."
    
    # -r is for single-end files. 
    # We use ${ACC}.fastq.gz because that's the output of fasterq-dump + pigz.
    salmon quant -i "$SALMON_INDEX" -l A \
         -r "${RAW_DIR}/${ACC}.fastq.gz" \
         --validateMappings \
         --threads $THREADS \
         -o "${OUT_DIR}/${ACC}_quant"
done

echo "[OK] Quantification complete. Results in $OUT_DIR"

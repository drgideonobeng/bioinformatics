#!/bin/bash
source config.sh

echo "=== Building Salmon Index ==="
echo "Using Transcriptome: $TRANSCRIPTOME_FA"
echo "Output Index: $SALMON_INDEX"

# The -i flag will create the directory automatically
salmon index \
    -t "$TRANSCRIPTOME_FA" \
    -i "$SALMON_INDEX" \
    --threads $THREADS

echo "----------------------------------------------------------------"
echo "[OK] Salmon index built. You are now ready for Quantification."
echo "----------------------------------------------------------------"

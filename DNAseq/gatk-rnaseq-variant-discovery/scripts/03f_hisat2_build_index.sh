#!/usr/bin/env bash

# Step 3f: Build HISAT2 genome index for GRCh38

source ./scripts/config.sh

HISAT2_INDEX_DIR="${REF_DIR}/hisat2_index_grch38"
mkdir -p "$HISAT2_INDEX_DIR"

# HISAT2 index base name (will create multiple .ht2 files)
HISAT2_INDEX_BASE="${HISAT2_INDEX_DIR}/GRCh38"

hisat2-build -p "$THREADS" "$FASTA" "$HISAT2_INDEX_BASE"

echo
echo "Index files:"
ls -lh "${HISAT2_INDEX_DIR}" | head

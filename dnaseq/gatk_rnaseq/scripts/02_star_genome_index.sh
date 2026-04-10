#!/usr/bin/env bash

# Build the STAR genome index (run from repo root)

# Load project settings and paths
source ./scripts/config.sh

# Make sure the index folder exists
mkdir -p "$STAR_INDEX"

# Build the index
"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX" \
  --genomeFastaFiles "$FASTA" \
  --sjdbGTFfile "$GTF" \
  --sjdbOverhang "$SJDB_OVERHANG"

# Show a few index files so you know it worked
ls -lh "$STAR_INDEX" | head

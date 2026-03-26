#!/usr/bin/env bash

# Step 3: Align RNA-seq FASTQs with STAR (single-pass)
# Output: coordinate-sorted BAM + logs

source ./scripts/config.sh

ALIGN_DIR="${RESULTS_DIR}/03_star_align_single_pass"
mkdir -p "$ALIGN_DIR"
OUT_PREFIX="${ALIGN_DIR}/SRR3192657."

# Clean old failed outputs (safe)
rm -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
rm -rf "${OUT_PREFIX}_STARtmp"

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" "$R2" \
  --outFileNamePrefix "$OUT_PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMunmapped Within

echo
echo "Main outputs:"
ls -lh "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" 2>/dev/null
ls -lh "${OUT_PREFIX}Log.final.out" 2>/dev/null

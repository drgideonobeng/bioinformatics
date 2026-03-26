#!/usr/bin/env bash

# Step 3c: STAR alignment (single-pass) using cleaned FASTQ headers
# Output: coordinate-sorted BAM + logs

source ./scripts/config.sh

ALIGN_DIR="${RESULTS_DIR}/03c_star_align_single_pass_clean_fastq"
mkdir -p "$ALIGN_DIR"
OUT_PREFIX="${ALIGN_DIR}/SRR3192657."

# Clean FASTQs made in Step 3b
R1_CLEAN="${ROOT_DIR}/data/intermediate/SRR3192657_1.clean.fastq"
R2_CLEAN="${ROOT_DIR}/data/intermediate/SRR3192657_2.clean.fastq"

# Clean old outputs (safe)
rm -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
rm -rf "${OUT_PREFIX}_STARtmp"

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1_CLEAN" "$R2_CLEAN" \
  --outFileNamePrefix "$OUT_PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif

echo
echo "Main outputs:"
ls -lh "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" 2>/dev/null
ls -lh "${OUT_PREFIX}Log.final.out" 2>/dev/null

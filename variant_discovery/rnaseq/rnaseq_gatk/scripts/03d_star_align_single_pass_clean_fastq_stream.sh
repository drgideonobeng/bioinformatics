#!/usr/bin/env bash

# Step 3d: STAR alignment (single-pass) using cleaned FASTQs,
#          read via a shell command (workaround for "0 input reads").

source ./scripts/config.sh

ALIGN_DIR="${RESULTS_DIR}/03d_star_align_single_pass_clean_fastq_stream"
mkdir -p "$ALIGN_DIR"
OUT_PREFIX="${ALIGN_DIR}/SRR3192657."

R1_CLEAN="${ROOT_DIR}/data/intermediate/SRR3192657_1.clean.fastq"
R2_CLEAN="${ROOT_DIR}/data/intermediate/SRR3192657_2.clean.fastq"

rm -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam"
rm -rf "${OUT_PREFIX}_STARtmp"

"$STAR_BIN" \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1_CLEAN" "$R2_CLEAN" \
  --readFilesCommand "/bin/bash" "-lc" \
  --outFileNamePrefix "$OUT_PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outBAMsortingBinsN 32 \
  --limitBAMsortRAM 12000000000

echo
echo "Check mapping summary:"
head -n 20 "${OUT_PREFIX}Log.final.out"

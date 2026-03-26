#!/usr/bin/env bash
set -euo pipefail

# Load configuration
source "$(dirname "$0")/config.sh"

# Choose input BAM:
# Direct input from the data folder (since we skipped align/markdup)
INPUT_BAM="${DATA_DIR}/${INPUT_BAM}"

echo "== Variant calling with HaplotypeCaller (GVCF mode) =="
echo "Input BAM: ${INPUT_BAM}"
echo "Output GVCF: ${GVCF}"

# Start the Timer
START_TIME=$(date +%s)

gatk HaplotypeCaller \
  --java-options "-Xmx10g" \
  -R "${REF_FASTA}" \
  -I "${INPUT_BAM}" \
  -O "${GVCF}" \
  -ERC GVCF \
  -L chr20 \
  --native-pair-hmm-threads ${THREADS}

# End the Timer
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo "HaplotypeCaller complete."
echo "Time taken: $(($DURATION / 60)) minutes and $(($DURATION % 60)) seconds."

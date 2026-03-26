#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "== Converting SAM -> BAM and sorting =="
echo "Input:  ${SAM}"
echo "Output: ${SORTED_BAM}"

# samtools view: convert SAM to BAM
# samtools sort: sort alignments by genomic coordinate
samtools view -bS "${SAM}" \
  | samtools sort -@ "${THREADS}" -o "${SORTED_BAM}"

# Create BAM index
samtools index "${SORTED_BAM}"

echo
echo "Sorted BAM created."
echo "Tip: You can delete the SAM to save space:"
echo "rm -f ${SAM}"

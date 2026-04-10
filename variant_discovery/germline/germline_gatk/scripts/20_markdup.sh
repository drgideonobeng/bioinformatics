#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "== Marking duplicates =="
echo "Input:  ${SORTED_BAM}"
echo "Output: ${MARKDUP_BAM}"
echo "Metrics: ${MARKDUP_METRICS}"

gatk MarkDuplicates \
  -I "${SORTED_BAM}" \
  -O "${MARKDUP_BAM}" \
  -M "${MARKDUP_METRICS}" \
  --CREATE_INDEX true

echo "MarkDuplicates complete."

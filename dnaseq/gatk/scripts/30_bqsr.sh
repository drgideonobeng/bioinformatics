#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "== BQSR: checking known-sites files =="
if [[ ! -f "${DBSNP_VCF}" ]]; then
  echo "WARNING: Missing DBSNP_VCF: ${DBSNP_VCF}"
  echo "BQSR will not run. Download known-sites or skip this step."
  exit 0
fi

if [[ ! -f "${INDELS_VCF}" ]]; then
  echo "WARNING: Missing INDELS_VCF: ${INDELS_VCF}"
  echo "BQSR will not run. Download known-sites or skip this step."
  exit 0
fi

echo
echo "== Building recalibration table =="
gatk BaseRecalibrator \
  -R "${REF_FASTA}" \
  -I "${MARKDUP_BAM}" \
  --known-sites "${DBSNP_VCF}" \
  --known-sites "${INDELS_VCF}" \
  -O "${RECAL_TABLE}"

echo
echo "== Applying BQSR =="
gatk ApplyBQSR \
  -R "${REF_FASTA}" \
  -I "${MARKDUP_BAM}" \
  --bqsr-recal-file "${RECAL_TABLE}" \
  -O "${RECAL_BAM}"

# Index recalibrated BAM
samtools index "${RECAL_BAM}"

echo "BQSR complete."

#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "== Genotyping GVCF -> VCF =="
echo "Input GVCF: ${GVCF}"
echo "Output VCF: ${RAW_VCF}"

gatk GenotypeGVCFs \
  -R "${REF_FASTA}" \
  -V "${GVCF}" \
  -O "${RAW_VCF}"

echo "GenotypeGVCFs complete."

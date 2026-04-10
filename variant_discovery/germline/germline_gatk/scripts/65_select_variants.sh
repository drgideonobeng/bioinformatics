#!/usr/bin/env bash 
set -euo pipefail 

# Load configuration source 
source "$(dirname "$0")/config.sh" 

echo "== selecting variants (PASS-only) =="
echo "Input VCF: ${FILTERED_VCF}" 
echo "Output VCF: ${PASS_VCF}"

gatk SelectVariants \
  -R "${REF_FASTA}" \
  -V "${FILTERED_VCF}" \
  -O "${PASS_VCF}"

echo
echo "== SelectVariants complete =="
echo "PASS count:"
zgrep -vc "^#" "${PASS_VCF}" || true

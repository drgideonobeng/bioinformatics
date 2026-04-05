#!/usr/bin/env bash
set -euo pipefail

# Load configuration
source "$(dirname "$0")/config.sh"

echo "== Hard-filtering variants (beginner practice defaults) =="
echo "Input VCF:  ${RAW_VCF}"
echo "Output VCF: ${FILTERED_VCF}"

# These filters label variants as suspicious; they do NOT delete them.
# PASS = not flagged by filters. FILTER field holds filter names if flagged.
gatk VariantFiltration \
  -R "${REF_FASTA}" \
  -V "${RAW_VCF}" \
  -O "${FILTERED_VCF}" \
  --filter-expression "QD < 2.0"  --filter-name "LowQD" \
  --filter-expression "FS > 60.0" --filter-name "HighFS" \
  --filter-expression "MQ < 30.0" --filter-name "LowMQ" \
  --filter-expression "DP < 3.0"  --filter-name "LowDP" \
  --genotype-filter-expression "GQ < 5.0" --genotype-filter-name "LowGQ" \
  --filter-expression "QUAL < 20.0" --filter-name "LowConfidence"

echo
echo "Filtering complete."
echo
echo "== Filter summary =="
# prints the PASS vs flagged breakdown
zgrep -v "^#" "${FILTERED_VCF}" | cut -f7 | sort | uniq -c || true

echo "Tip: view first non-header lines:"
echo "zgrep -v '^#' ${FILTERED_VCF} | head"

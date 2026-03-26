#!/usr/bin/env bash
set -euo pipefail

# Runs the pipeline in order.
# If BQSR can't run (missing known-sites), it will skip that step safely.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"${SCRIPT_DIR}/00_check_setup.sh"
"${SCRIPT_DIR}/05_index_reference.sh"
"${SCRIPT_DIR}/10_align.sh"
"${SCRIPT_DIR}/15_sort_bam.sh"
"${SCRIPT_DIR}/20_markdup.sh"
"${SCRIPT_DIR}/30_bqsr.sh"
"${SCRIPT_DIR}/40_haplotypecaller.sh"
"${SCRIPT_DIR}/45_genotype_gvcf.sh"
"${SCRIPT_DIR}/60_filter_variants.sh"

echo
echo "All done. Final output:"
echo "results/*.filtered.vcf.gz"

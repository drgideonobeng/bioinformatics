#!/usr/bin/env bash
# set -eou pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 4: Statistical Testing (MAGeCK Test)"
echo "======================================================"

# The normalized counts file generated in Step 3
COUNTS_MATRIX="${COUNT_DIR}/${PROJECT_NAME}.count.txt"

echo "=> Running MAGeCK test (Treatment vs. Control)..."
mageck test \
    -k "${COUNTS_MATRIX}" \
    -t "${TREAT_SAMPLES}" \
    -c "${CTRL_SAMPLES}" \
    -n "${TEST_DIR}/${PROJECT_NAME}_drug_vs_veh" \
    --norm-method median \
    --pdf-report

echo "=> Top 10 Negatively Selected Genes (Drug Sensitizers):"
head -n 11 "${TEST_DIR}/${PROJECT_NAME}_drug_vs_veh.gene_summary.txt" | awk '{print $1, $2, $3, $8, $9}' | column -t

echo "======================================================"
echo " Pipeline Complete. Check the ${TEST_DIR} folder."
echo "======================================================"

#!/usr/bin/env bash
# set -eou pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 3: Quantifying sgRNAs (MAGeCK Count)"
echo "======================================================"

# 1. Prepare a SPACE-separated list of all trimmed FASTQ files
# MAGeCK needs spaces between different files to count them as separate samples
FASTQ_LIST=$(ls "${TRIM_DIR}"/*_trimmed.fastq.gz | xargs)

# 2. Prepare a COMMA-separated list of sample names 
SAMPLE_NAMES=$(ls "${TRIM_DIR}"/*_trimmed.fastq.gz | xargs -n 1 basename | sed 's/_trimmed.fastq.gz//g' | paste -sd "," -)

echo "=> Running MAGeCK count..."
# Note: We use ${FASTQ_LIST} without quotes so the shell expands the spaces into separate arguments
mageck count \
    -l "${LIBRARY_FILE}" \
    --fastq ${FASTQ_LIST} \
    --sample-label "${SAMPLE_NAMES}" \
    -n "${COUNT_DIR}/${PROJECT_NAME}" \
    --trim-5 0 \
    --pdf-report

echo "=> Count summary (First 10 lines):"
head -n 10 "${COUNT_DIR}/${PROJECT_NAME}.countsummary.txt" | column -t

echo "=> Quantification complete!"

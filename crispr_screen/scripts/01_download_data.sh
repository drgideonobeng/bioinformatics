#!/usr/bin/env bash
# set -eou pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 1: Downloading Raw Data & Library"
echo "======================================================"

echo "=> Downloading library reference..."
if [ ! -f "${LIBRARY_FILE}" ]; then
    wget -q -O "${LIBRARY_FILE}" "${LIBRARY_URL}"
fi

echo "=> Downloading FASTQ files..."
declare -A SAMPLES=(
    ["T0-Control"]="${T0_CTRL_URL}"
    ["T8-Vehicle"]="${T8_VEH_URL}"
    ["T8-APR-246"]="${T8_DRUG_URL}"
)

for SAMP in "${!SAMPLES[@]}"; do
    FILE_PATH="${RAW_DIR}/${SAMP}.fastq.gz"
    if [ ! -f "${FILE_PATH}" ]; then
        echo "   -> Fetching ${SAMP}..."
        wget -q -O "${FILE_PATH}" "${SAMPLES[$SAMP]}"
    else
        echo "   -> ${SAMP} already exists. Skipping."
    fi
done

echo "=> Download complete!"

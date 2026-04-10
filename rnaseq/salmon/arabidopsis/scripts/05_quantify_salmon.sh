#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 5: Quantifying Clean Expression with Salmon"
echo "======================================================"

# Loop through the exact same SRR_LIST used in the download and trim steps
for SRR in $(cat "${SRR_LIST}"); do
    echo "=> Processing sample: ${SRR}..."
    salmon quant -i "${INDEX_DIR}" \
                 -l "${LIB_TYPE}" \
                 -1 "${TRIM_DIR}/${SRR}/${SRR}_1.trim.fastq.gz" \
                 -2 "${TRIM_DIR}/${SRR}/${SRR}_2.trim.fastq.gz" \
                 -p "${THREADS}" \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 -o "${RESULTS_DIR}/${SRR}"
    echo "   Finished ${SRR}."
done

echo "===== Showing the first 10 lines of each quantification file ======"

for SRR in $(cat "${SRR_LIST}"); do
    echo "--- ${SRR} ---"
    head "${RESULTS_DIR}/${SRR}/quant.sf" | column -t
    echo "" 
done

echo "======================================================"
echo " Step 5: Quantification completed"
echo "======================================================"

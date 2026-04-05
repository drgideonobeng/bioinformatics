#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

# Define a new directory for the clean data
export TRIM_DIR="${ROOT_DIR}/data_trimmed"
mkdir -p "${TRIM_DIR}"

echo "======================================================"
echo " Step 4: Trimming Adapters & Quality Filtering (fastp)"
echo "======================================================"

# Loop through each SRR number in your text file
for SRR in $(cat "${SRR_LIST}"); do
    echo "------------------------------------------------"
    echo " Trimming sample: ${SRR}"
    
    # Define input paths (your raw data)
    IN1="${DATA_DIR}/${SRR}/${SRR}_1.fastq.gz"
    IN2="${DATA_DIR}/${SRR}/${SRR}_2.fastq.gz"
    
    # Create a subfolder for the trimmed sample
    OUT_SAMP_DIR="${TRIM_DIR}/${SRR}"
    mkdir -p "${OUT_SAMP_DIR}"
    
    # Define output paths (your clean data)
    OUT1="${OUT_SAMP_DIR}/${SRR}_1.trim.fastq.gz"
    OUT2="${OUT_SAMP_DIR}/${SRR}_2.trim.fastq.gz"
    
    # Skip if already trimmed
    if [ -f "${OUT1}" ] && [ -f "${OUT2}" ]; then
        echo "   -> Trimmed files already exist. Skipping."
        continue
    fi

    # Run fastp
    # --detect_adapter_for_pe: Automatically finds paired-end Illumina adapters
    # --thread: Uses the threads from your config
    fastp -i "${IN1}" -I "${IN2}" \
          -o "${OUT1}" -O "${OUT2}" \
          --detect_adapter_for_pe \
          --thread ${THREADS} \
          --html "${OUT_SAMP_DIR}/${SRR}_fastp.html" \
          --json "${OUT_SAMP_DIR}/${SRR}_fastp.json" 2> "${OUT_SAMP_DIR}/${SRR}_fastp.log"
          
    echo "   -> ${SRR} complete!"
done

echo "======================================================"
echo "=> All data perfectly trimmed and saved in: ${TRIM_DIR}"
echo "======================================================"

#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

# Define and create a specific directory for QC reports
QC_DIR="${OUT_DIR}/fastqc"
mkdir -p "${QC_DIR}"

echo "======================================================"
echo " Step 3: Quality Control (FastQC & MultiQC)"
echo "======================================================"

echo "=> Running FastQC on all downloaded samples..."
# Run FastQC on all fastq.gz files in the data directory simultaneously
# It automatically uses the $THREADS variable from config.sh
fastqc -q -o "${QC_DIR}" -t "${THREADS}" "${DATA_DIR}"/*/*.fastq.gz

echo "=> Aggregating QC reports with MultiQC..."
# Run MultiQC to combine all FastQC reports into one interactive HTML file
multiqc -q -f -o "${QC_DIR}" "${QC_DIR}"

echo "======================================================"
echo "=> QC complete! Open ${QC_DIR}/multiqc_report.html to view."
echo "======================================================"

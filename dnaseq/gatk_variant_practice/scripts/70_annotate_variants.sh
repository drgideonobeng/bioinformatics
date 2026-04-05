#!/bin/bash
set -e

source scripts/config.sh

# ==========================================
# CONFIGURATION
# ==========================================
# NOTE: You must have the Funcotator data sources downloaded!
# If you don't have them, download the "small" version from GATK's FTP.
# Update this path to where you extracted the data source:
FUNC_DATA_SOURCE="${ROOT_DIR}/funcotator_dataSources.v1.7.20200521s"

ANNOTATED_VCF="${RESULTS_DIR}/${SAMPLE}.annotated.vcf"
ANNOTATED_TABLE="${RESULTS_DIR}/${SAMPLE}.annotated.table.txt"

echo "=========================================="
echo "Running Funcotator (Annotation)"
echo "Input VCF: ${PASS_VCF}"
echo "Output VCF: ${ANNOTATED_VCF}"
echo "Using Data Source: ${FUNC_DATA_SOURCE}"
echo "=========================================="

# Check if data source exists before running
if [ ! -d "${FUNC_DATA_SOURCE}" ]; then
    echo "ERROR: Funcotator data source not found at ${FUNC_DATA_SOURCE}"
    echo "Please download it or update the path in this script."
    exit 1
fi

gatk Funcotator \
    -R "${REF_FASTA}" \
    -V "${{ANNOTATED_VCF}" \
    --output-file-format VCF \
    --data-sources-path "${FUNC_DATA_SOURCE}" \
    --ref-version hg38 \
    --remove-filtered-variants

# Create a readable table (optional but nice)
# This extracts just the gene name and protein change for easier reading
echo "Creating summary table..."
gatk VariantsToTable \
    -R "${REF_FASTA}" \
    -V "${ANNOTATED_VCF}" \
    -F CHROM -F POS -F REF -F ALT -F FUNCOTATION \
    -O "${ANNOTATED_TABLE}"

echo "=========================================="
echo "Annotation complete!"
echo "Check the readable table here: ${ANNOTATED_TABLE}"
echo "=========================================="

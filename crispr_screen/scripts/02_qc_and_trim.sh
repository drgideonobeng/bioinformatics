#!/usr/bin/env bash
# set -eou pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 2: Quality Control & Adapter Trimming"
echo "======================================================"

echo "=> Running FastQC on raw reads..."
fastqc -q -o "${QC_DIR}" -t "${THREADS}" "${RAW_DIR}"/*.fastq.gz

echo "=> Trimming 5' adapters with Cutadapt..."
for FASTQ in "${RAW_DIR}"/*.fastq.gz; do
    BASENAME=$(basename "${FASTQ}" .fastq.gz)
    TRIMMED="${TRIM_DIR}/${BASENAME}_trimmed.fastq.gz"
    
    if [ ! -f "${TRIMMED}" ]; then
        echo "   -> Trimming ${BASENAME}..."
        # -g flags the 5' adapter. --discard-untrimmed removes reads without adapters
        cutadapt -g "${ADAPTER_5P}" \
                 -o "${TRIMMED}" \
                 -j "${THREADS}" \
                 --minimum-length 18 \
                 --discard-untrimmed \
                 "${FASTQ}" > "${QC_DIR}/${BASENAME}_cutadapt.log"
    fi
done

echo "=> Aggregating QC metrics with MultiQC..."
multiqc -q -f -o "${QC_DIR}" "${QC_DIR}"

echo "=> Trimming & QC complete!"

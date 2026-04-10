#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 5: Aligning Reads with STAR"
echo "======================================================"

# Loop through each sample in your list
for SRR in $(cat "${SRR_LIST}"); do
    echo "--------------------------------------------------"
    echo "=> Aligning sample: ${SRR}..."

    OUT_PREFIX="${ALIGN_DIR}/${SRR}_"

    # Check if BAM already exists to prevent redundant runs
    if [ -f "${OUT_PREFIX}Aligned.sortedByCoord.out.bam" ]; then
        echo "   => BAM file already exists. Skipping alignment."
        continue
        fi

# Run STAR alignment
    STAR --runThreadN "${THREADS}" \
         --genomeDir "${STAR_INDEX_DIR}" \
         --readFilesIn "${TRIM_DIR}/${SRR}/${SRR}_1.trim.fastq.gz" "${TRIM_DIR}/${SRR}/${SRR}_2.trim.fastq.gz" \
         --readFilesCommand "gunzip -c" \
         --outFileNamePrefix "${OUT_PREFIX}" \
         --outSAMtype BAM SortedByCoordinate \
         --outBAMsortingThreadN "${THREADS}" \
         --limitBAMsortRAM 4000000000 # Allocates RAM per thread for sorting

    echo "   Finished aligning ${SRR}."
done

echo "======================================================"
echo "=> All alignments complete! BAM files saved to: ${ALIGN_DIR}"
echo "======================================================"

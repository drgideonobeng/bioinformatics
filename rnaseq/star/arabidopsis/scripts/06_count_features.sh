#!/usr/bin/env bash

#Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "====================================================="
echo " Step 6: Quantifying Genes with featureCounts"
echo "====================================================="

FINAL_MATRIX="${COUNTS_DIR}/featureCounts_matrix.txt"

if [ -f "${FINAL_MATRIX}" ]; then
    echo "=> Count matrix already exists at ${FINAL_MATRIX}. Skipping."
else
    echo "=> Running featureCounts on all BAM files..."
    
    # -T: threads
    # -p: specifies that the data is paired-end
    # -a: path to the GTF annotation file
    # -o: output file path
    # We use a wildcard (*) to pass all sorted BAM files at once
    featureCounts -T "${THREADS}" \
                  -p \
                  -a "${ANNOTATION_GTF}" \
                  -o "${FINAL_MATRIX}" \
                  "${ALIGN_DIR}"/*Aligned.sortedByCoord.out.bam
                  
    echo "=> Feature counting complete!"
fi

echo ""
echo "===== Preview of the Final Count Matrix ======"
# The output has a comment line (1), header line (2), and then data.
# We cut to show column 1 (GeneID) and columns 7+ (which contain the actual sample counts)
head -n 10 "${FINAL_MATRIX}" | cut -f1,7- | column -t

echo "======================================================"
echo "=> Matrix ready for DESeq2 at: ${FINAL_MATRIX}"
echo "======================================================"

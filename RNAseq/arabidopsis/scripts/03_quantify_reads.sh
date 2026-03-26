#!/usr/bin/env bash

# Load the master configuration file
source "$(dirname "$0")/../config.sh"

echo "======================================================"
echo " Step 3: Quantifying Expression with Salmon"
echo "======================================================"

for SAMP in "${SAMPLES[@]}"; do
    echo "=> Processing sample: ${SAMP}..."
    
    # Run Salmon in quantification mode
    salmon quant -i "${INDEX_DIR}" \
                 -l "${LIB_TYPE}" \
                 -1 "${DATA_DIR}/${SAMP}/${SAMP}_1.fastq.gz" \
                 -2 "${DATA_DIR}/${SAMP}/${SAMP}_2.fastq.gz" \
                 -p "${THREADS}" \
                 --validateMappings \
                 -o "${RESULTS_DIR}/${SAMP}_quant"
                 
    echo "   Finished ${SAMP}."
done

echo "======================================================"
echo "=> Quantification complete!"

echo "=> Showing the first 10 lines of each sample"
head -n 10 results/quants/DRR016125_quant/quant.sf

echo "=> Results are saved in: ${RESULTS_DIR}"
echo "=====Step 03_Quantification is complete"

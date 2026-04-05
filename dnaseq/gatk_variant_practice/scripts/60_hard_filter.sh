#!/bin/bash
set -e

source scripts/config.sh

# 5. Merge SNPs and Indels back together
echo "Step 5: Merging filtered SNPs and Indels..."
gatk MergeVcfs \
    -I "${RESULTS_DIR}/${SAMPLE}.snps.filtered.vcf.gz" \
    -I "${RESULTS_DIR}/${SAMPLE}.indels.filtered.vcf.gz" \
    -O "${FILTERED_VCF}"

# 6. Create the Final "PASS" Only File
echo "Step 6: Keeping only PASS variants..."
gatk SelectVariants \
    -R "${REF_FASTA}" \
    -V "${FILTERED_VCF}" \
    --exclude-filtered \
    -O "${PASS_VCF}"

echo "=========================================="
echo "Filtering complete!"
echo "Final 'Good' Variants: ${PASS_VCF}"
echo "=========================================="

# 7. Summary Check
echo "Filter Distribution:"
zgrep -v "^#" "${FILTERED_VCF}" | cut -f7 | sort | uniq -c

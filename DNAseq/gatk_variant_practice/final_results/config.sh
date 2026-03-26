#!/usr/bin/env bash
set -euo pipefail

# ========== BASIC PROJECT PATHS ==========
# ROOT_DIR is the repository root (one level above scripts/)
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA_DIR="${ROOT_DIR}/data"
REF_DIR="${ROOT_DIR}/ref"
KNOWN_DIR="${ROOT_DIR}/known-sites"
RESULTS_DIR="${ROOT_DIR}/results_full_chr20"

mkdir -p "${RESULTS_DIR}"

# ========== SAMPLE INPUTS ==========
# Update these filenames to match your actual FASTQs
SAMPLE="NA12878_chr20_subset"
R1="${DATA_DIR}/H06HDADXX130110.1.ATCACGAT.20k_reads_1.fastq"
R2="${DATA_DIR}/H06HDADXX130110.1.ATCACGAT.20k_reads_2.fastq"
INPUT_BAM="NA12878_24RG_small.hg38.bam"

# ========== REFERENCE ==========
# Update to your reference fasta name
REF_FASTA="${REF_DIR}/Homo_sapiens_assembly38.fasta"

# ========== KNOWN SITES (for BQSR) ==========
# Optional, but recommended if you have them
DBSNP_VCF="${KNOWN_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"
INDELS_VCF="${KNOWN_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ========== READ GROUP (required by GATK) ==========
# Keep simple for practice
RGID="${SAMPLE}_chr20"
RGLB="lib1"
RGPL="ILLUMINA"
RGPU="unit1"
RGSM="${SAMPLE}"

READ_GROUP="@RG\tID:${RGID}\tSM:${RGSM}\tPL:${RGPL}\tLB:${RGLB}\tPU:${RGPU}"

# ========== THREADS ==========
# Adjust if you want
THREADS=4

# ========== OUTPUT FILES ==========
SAM="${RESULTS_DIR}/${SAMPLE}.sam"
SORTED_BAM="${RESULTS_DIR}/${SAMPLE}.sorted.bam"
MARKDUP_BAM="${RESULTS_DIR}/${SAMPLE}.markdup.bam"
MARKDUP_METRICS="${RESULTS_DIR}/${SAMPLE}.markdup.metrics.txt"

RECAL_TABLE="${RESULTS_DIR}/${SAMPLE}.recal.table"
RECAL_BAM="${RESULTS_DIR}/${SAMPLE}.recal.bam"

export GVCF="${RESULTS_DIR}/${SAMPLE}.g.vcf.gz"
export RAW_VCF="${RESULTS_DIR}/${SAMPLE}.raw.vcf.gz"
export FILTERED_VCF="${RESULTS_DIR}/${SAMPLE}.filtered.vcf.gz"
export PASS_VCF="${RESULTS_DIR}/${SAMPLE}.pass.vcf.gz"


#!/usr/bin/env bash

# ========== BASIC PROJECT PATHS ==========
# ROOT_DIR is the repository root (one level above scripts/)
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ========== SETTINGS ==========
THREADS=8
SJDB_OVERHANG=149   # usually read_length - 1 (150bp reads -> 149)

# ========== PROJECT FOLDERS ==========
SCRIPTS_DIR="${ROOT_DIR}/scripts"
RAW_DIR="${ROOT_DIR}/data/raw"
REF_DIR="${ROOT_DIR}/data/ref"
RESULTS_DIR="${ROOT_DIR}/results"

# ========== INPUT FASTQs (paired-end) ==========
R1="${RAW_DIR}/SRR3192657_1.fastq"
R2="${RAW_DIR}/SRR3192657_2.fastq"

# ========== REFERENCE FILES ==========
FASTA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${REF_DIR}/gencode.v47.annotation.gtf"

# ========== STAR INDEX ==========
STAR_INDEX="${REF_DIR}/star_index_gencode_v47"

# ========== TOOLS ==========
# Use the STAR aligner installed inside your active conda env
STAR_BIN="${CONDA_PREFIX}/bin/STAR"

#!/usr/bin/env bash

# This script downloads one RNA-seq dataset (SRR3192657) from NCBI
# and converts it into two paired-end FASTQ files.

# 1) Download the dataset as an .sra file into data/raw/
prefetch SRR3192657 --output-directory data/raw

# 2) Convert the .sra file into paired FASTQ files
#    --split-files makes two files: _1 and _2
#    --outdir chooses where the FASTQs go
fasterq-dump data/raw/SRR3192657/SRR3192657.sra --split-files --outdir data/raw

# 3) Compress the FASTQ files to save space (STAR can read .gz files)
gzip -f data/raw/SRR3192657_1.fastq
gzip -f data/raw/SRR3192657_2.fastq

# 4) Show the final files
ls -lh data/raw/SRR3192657_*.fastq.gz

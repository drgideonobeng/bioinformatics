#!/usr/bin/env bash

# Step 3b: Fix FASTQ headers to be STAR-friendly for paired-end data
# Input:  data/raw/SRR3192657_1.fastq and _2.fastq
# Output: data/intermediate/SRR3192657_1.clean.fastq and _2.clean.fastq

source ./scripts/config.sh

mkdir -p "${ROOT_DIR}/data/intermediate"

R1_IN="$R1"
R2_IN="$R2"
R1_OUT="${ROOT_DIR}/data/intermediate/SRR3192657_1.clean.fastq"
R2_OUT="${ROOT_DIR}/data/intermediate/SRR3192657_2.clean.fastq"

echo "Writing:"
echo "  $R1_OUT"
echo "  $R2_OUT"

# For every FASTQ record:
# - line 1 is header: keep only first token and add /1 or /2
# - lines 2,3,4 are unchanged
awk 'NR%4==1{split($1,a,"@"); print "@"a[2]"/1"; next} {print}' "$R1_IN" > "$R1_OUT"
awk 'NR%4==1{split($1,a,"@"); print "@"a[2]"/2"; next} {print}' "$R2_IN" > "$R2_OUT"

echo
echo "Example headers after cleaning:"
head -n 1 "$R1_OUT"
head -n 1 "$R2_OUT"

echo
echo "Output sizes:"
ls -lh "$R1_OUT" "$R2_OUT"

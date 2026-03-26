#!/usr/bin/env bash
set -euo pipefail
# Load configuration
source "$(dirname "$0")/config.sh"

echo "== Indexing reference for samtools (FAI) =="
if [[ ! -f "${REF_FASTA}.fai" ]]; then
  samtools faidx "${REF_FASTA}"
else
  echo "Skipping: ${REF_FASTA}.fai already exists"
fi

echo
echo "== Creating GATK sequence dictionary (.dict) =="
DICT="${REF_FASTA%.*}.dict"
if [[ ! -f "${DICT}" ]]; then
  gatk CreateSequenceDictionary -R "${REF_FASTA}" -O "${DICT}"
else
  echo "Skipping: ${DICT} already exists"
fi

echo
echo "== Indexing reference for bwa =="
# BWA creates several files with extensions like .bwt, .pac, .ann, .amb, .sa
if [[ ! -f "${REF_FASTA}.bwt" ]]; then
  bwa index "${REF_FASTA}"
else
  echo "Skipping: bwa index files already exist"
fi

echo
echo "Reference indexing done."

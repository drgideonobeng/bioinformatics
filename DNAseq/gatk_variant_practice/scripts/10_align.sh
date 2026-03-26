#!/usr/bin/env bash
set -euo pipefail
# Load configuration
source "$(dirname "$0")/config.sh"

echo "== Aligning reads with bwa mem =="
echo "Output: ${SAM}"
echo "Read group: ${READ_GROUP}"

bwa mem -M -t "${THREADS}" -R "${READ_GROUP}" \
  "${REF_FASTA}" \
  "${R1}" "${R2}" \
  > "${SAM}"

echo "Alignment complete."

#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/config.sh"

echo "== Checking tools =="
command -v gatk >/dev/null && echo "OK: gatk found" || { echo "ERROR: gatk not found"; exit 1; }
command -v bwa  >/dev/null && echo "OK: bwa found"  || { echo "ERROR: bwa not found";  exit 1; }
command -v samtools >/dev/null && echo "OK: samtools found" || { echo "ERROR: samtools not found"; exit 1; }

echo
echo "== Checking files =="
[[ -f "${R1}" ]] && echo "OK: R1 exists: ${R1}" || { echo "ERROR: missing R1: ${R1}"; exit 1; }
[[ -f "${R2}" ]] && echo "OK: R2 exists: ${R2}" || { echo "ERROR: missing R2: ${R2}"; exit 1; }
[[ -f "${REF_FASTA}" ]] && echo "OK: REF exists: ${REF_FASTA}" || { echo "ERROR: missing REF: ${REF_FASTA}"; exit 1; }

echo
echo "== Folders =="
echo "DATA_DIR:    ${DATA_DIR}"
echo "REF_DIR:     ${REF_DIR}"
echo "KNOWN_DIR:   ${KNOWN_DIR}"
echo "RESULTS_DIR: ${RESULTS_DIR}"

echo
echo "Setup looks good."

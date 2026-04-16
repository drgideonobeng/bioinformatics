#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 1: Downloading References & Building Decoys"
echo "======================================================"

# --- 1. DOWNLOAD TRANSCRIPTOME ---
if [ ! -f "${TRANSCRIPTOME_FA}" ]; then
    echo "=> Downloading transcriptome..."
    curl -L -o "${TRANSCRIPTOME_FA}" "${TRANSCRIPTOME_URL}"
fi

# --- 2. DOWNLOAD GENOME (For Decoys) ---
if [ ! -f "${GENOME_FA}" ]; then
    echo "=> Downloading full genome for decoy indexing..."
    curl -L -o "${GENOME_FA}" "${GENOME_URL}"
fi

# --- 3. CREATE DECOY FILE & GENTROME ---
DECOYS="${REF_DIR}/decoys.txt"
GENTROME="${REF_DIR}/gentrome.fa.gz"

if [ ! -f "${GENTROME}" ]; then
    echo "=> Extracting chromosome names for decoys.txt..."
    # Unzip the genome on the fly, find headers (>), grab the first word, and remove the '>'
    gunzip -c "${GENOME_FA}" | grep "^>" | cut -d " " -f 1 | tr -d ">" > "${DECOYS}"

    echo "=> Concatenating transcriptome and genome into gentrome.fa.gz..."
    cat "${TRANSCRIPTOME_FA}" "${GENOME_FA}" > "${GENTROME}"
fi

echo "======================================================"
echo " Building the Decoy-Aware Salmon Index"
echo "======================================================"

# Check if the index is already built by looking for a core Salmon file
if [ -f "${INDEX_DIR}/seq.bin" ]; then
    echo "=> Salmon index already exists at ${INDEX_DIR}. Skipping..."
else
    echo "=> Generating Salmon Index (This may take a bit more RAM!)..."
    salmon index -t "${GENTROME}" -d "${DECOYS}" -i "${INDEX_DIR}" -k ${KMER_SIZE}
    echo "=> Decoy-aware index successfully built at: ${INDEX_DIR}"
fi

echo "======================================================"
echo " Step 01 : Downloading Ref and Indexing Complete"
echo "======================================================"

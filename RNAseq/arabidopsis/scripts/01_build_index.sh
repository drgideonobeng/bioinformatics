#!/usr/bin/env bash

# Load the master configuration file
source "$(dirname "$0")/../config.sh"

echo "======================================================"
echo " Step 1: Building the Salmon Index"
echo "======================================================"

echo "=> Downloading Transcriptome from Ensembl..."
curl -L -o "${TRANSCRIPTOME_FA}" "${TRANSCRIPTOME_URL}"

echo "=> Generating Salmon Index..."
# -t: transcripts, -i: output index, -k: k-mer size
salmon index -t "${TRANSCRIPTOME_FA}" -i "${INDEX_DIR}" -k ${KMER_SIZE}

echo "=> Index successfully built at: ${INDEX_DIR}"

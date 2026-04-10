#!/usr/bin/env bash

# Get the absolute path to the directory this script lives in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the config file located one folder up
source "${SCRIPT_DIR}/../config.sh"

echo "======================================================"
echo " Step 1: Downloading References & Building STAR Index"
echo "======================================================"

# --- 1. DOWNLOAD & UNZIP GENOME ---
if [ ! -f "${GENOME_FA}" ]; then
    echo "=> Downloading full genome..."
    curl -L -o "${GENOME_FA}.gz" "${GENOME_URL}"
    
    echo "=> Unzipping genome..."
    gunzip -f "${GENOME_FA}.gz"
fi

# --- 2. DOWNLOAD & UNZIP ANNOTATION (GTF) ---
if [ ! -f "${ANNOTATION_GTF}" ]; then
    echo "=> Downloading GTF annotation..."
    curl -L -o "${ANNOTATION_GTF}.gz" "${ANNOTATION_URL}"
    
    echo "=> Unzipping GTF..."
    gunzip -f "${ANNOTATION_GTF}.gz"
fi

# --- 3. BUILD STAR INDEX  ---
if [ -f "${STAR_INDEX_DIR}/Genome" ]; then
    echo "=> STAR index already exists at ${STAR_INDEX_DIR}. Skipping ..."
else 
echo "=> Generating STAR Index (This requires 30GB+ of RAM for large genomes, but ~3GB for Arabidopsis)..."
    STAR --runThreadN "${THREADS}" \
         --runMode genomeGenerate \
         --genomeDir "${STAR_INDEX_DIR}" \
         --genomeFastaFiles "${GENOME_FA}" \
         --sjdbGTFfile "${ANNOTATION_GTF}" \
         --sjdbOverhang "${OVERHANG}" \
         --genomeSAindexNbases 12

    echo "=> STAR index successfully built"
fi

#!/bin/bash
source config.sh

echo "=== Downloading Yeast Transcriptome (cDNA) ==="

# Download the cDNA FASTA
curl -L -o "${TRANSCRIPTOME_FA}.gz" ftp://ftp.ensembl.org/pub/release-104/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz

# Download the GTF (we will need this for R later!)
curl -L -o "${GENOME_GTF}.gz" ftp://ftp.ensembl.org/pub/release-104/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.104.gtf.gz

echo "=== Decompressing Files ==="
gunzip -f "${TRANSCRIPTOME_FA}.gz"
gunzip -f "${GENOME_GTF}.gz"

echo "----------------------------------------------------------------"
echo "[OK] Reference files are now physically in $REF_DIR"
ls -lh "$REF_DIR"
echo "----------------------------------------------------------------"

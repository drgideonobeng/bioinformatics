# CRISPR Screen Analysis Pipeline

This repository contains an automated, industry-standard pipeline for analyzing pooled CRISPR-Cas9 knockout screens. The workflow processes raw FASTQ files to identify essential genes driving phenotypes (e.g., drug resistance or sensitivity) using the MAGeCK algorithm.

## Workflow Overview
1. **Data Acquisition:** Downloads raw FASTQ sequences and the sgRNA library reference.
2. **Quality Control & Trimming:** Uses `FastQC`, `MultiQC`, and `Cutadapt` to remove vector adapter sequences.
3. **Quantification:** Maps reads to the guide library using `mageck count`.
4. **Statistical Testing:** Identifies enriched/depleted genes using `mageck test` (Robust Rank Aggregation).

## Quick Start
1. Install Conda/Mamba and create the environment:
   `conda env create -f env.yml`
   `conda activate crispr_pipeline`
2. Adjust parameters in `config.sh` if running a different dataset.
3. Run the scripts sequentially from the project root:
   `bash scripts/01_download_data.sh`
   `bash scripts/02_qc_and_trim.sh`
   `bash scripts/03_mageck_count.sh`
   `bash scripts/04_mageck_test.sh`

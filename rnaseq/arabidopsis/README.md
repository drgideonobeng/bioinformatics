# Arabidopsis RNA-Seq Quantification: A Salmon Pseudoalignment Pipeline

This repository contains a reproducible, production-level bioinformatics pipeline for processing bulk RNA-seq data. It utilizes **Salmon** for fast, bias-aware, transcript-level quantification using pseudoalignment.

The data used in this project originates from an *Arabidopsis thaliana* RNA-seq experiment ([PRJDB2508](https://combine-lab.github.io/salmon/getting_started/)). 

## Why Salmon? (Pseudoalignment vs. Traditional Alignment)
Traditional RNA-seq pipelines (e.g., HISAT2 + featureCounts) map reads against the entire genome, which is computationally expensive and generates massive intermediate BAM files. 

This pipeline leverages **Salmon**, which maps reads directly to the **transcriptome** (cDNA). This pseudoalignment approach offers several advantages:
* **Speed:** Quantifies samples in minutes rather than hours.
* **Accuracy:** Automatically corrects for GC-content and sequence-specific biases.
* **Efficiency:** Eliminates the need for massive intermediate disk storage.

## Project Architecture
The project is built with modularity and reproducibility in mind, utilizing a centralized `config.sh` to manage paths and parameters, preventing hard-coded fragility.

```text
arabidopsis-salmon-rnaseq/
├── README.md              # Project documentation
├── environment.yml        # Conda environment definition
├── config.sh              # Centralized global variables and paths
├── scripts/
│   ├── 01_build_index.sh  # Downloads cDNA and builds the Salmon index
│   ├── 02_download.sh     # Robust paired-end read acquisition via SRA
│   └── 03_quantify.sh     # Core Salmon pseudoalignment and quantification
├── data/                  # Raw FASTQ files (git-ignored)
├── reference/             # Transcriptome FASTA and Salmon Index (git-ignored)
└── results/
    └── quants/            # Final transcript counts (.sf files) ready for R

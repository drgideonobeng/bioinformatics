# Bioinformatics Portfolio

Welcome to my bioinformatics portfolio! This repository contains end-to-end computational biology pipelines (automated), data analysis scripts, and visualization workflows. 

These projects demonstrate my ability to process raw high-throughput sequencing data, perform rigorous quality control, execute statistical analyses, and generate publication-ready visualizations.

---

## ️ Tech Stack & Tools

* **Languages:** R, Python, Bash / Shell Scripting
* **Transcriptomics (Bulk RNA-seq):** DESeq2, Salmon, STAR, HISAT2, clusterProfiler
* **Single-Cell (scRNA-seq):** Seurat (v4/v5), dimensionality reduction (PCA, UMAP), marker identification
* **Genomics (DNA-seq):** GATK (Genome Analysis Toolkit), Variant Calling (BWA, HaplotypeCaller)
* **Quality Control:** FastQC, MultiQC
* **Environment Management:** Conda (`environment.yml`)

---

## Repository Structure & Projects

### 1. Single-Cell RNA-Seq Analysis (`scRNAseq/`)
End-to-end processing of 10x Genomics single-cell data, including QC, normalization, clustering, and cell type annotation using **Seurat**.
* **`pbmc-3k/` & `pbmc-10k/`:** Standard clustering and annotation workflows for peripheral blood mononuclear cells.
* **`mouse-brain/`:** scRNA-seq analysis of complex tissue types.
* **`pbmc-multisample/`:** Advanced pipeline integrating multiple samples to evaluate differential expression between treated and untreated conditions (e.g., B-cell response pathways).

### 2. Bulk RNA-Seq Transcriptomics (`RNAseq/`)
Pipelines for quantifying gene expression, performing differential expression analysis (DEA), and pathway enrichment.
* **`airway-human-rnaseq/`:** Differential expression analysis using **DESeq2** with custom GSEA and network plotting.
* **`arabidopsis/` & `yeast-rnaseq/`:** Handling non-human transcriptomes, pseudo-alignment with **Salmon**, and MultiQC integration.
* **`rnaseq-colon-cancer-analysis/`:** In-depth cancer transcriptomics, featuring PCA biplots, survival analysis, and driver gene identification.

### 3. DNA-Seq Variant Discovery (`DNAseq/`)
Implementation of the Broad Institute's Best Practices for variant calling.
* **`gatk_variant_practice/`:** Full pipeline for marking duplicates, Base Quality Score Recalibration (BQSR), and calling SNPs/Indels using **GATK HaplotypeCaller**.
* **`gatk-rnaseq-variant-discovery/`:** Calling variants directly from RNA-seq reads, including STAR 2-pass alignment and read formatting.

### 4. Data Wrangling & Evaluation (`target-eval/`)
* Python-based scripts for cleaning messy expression matrices, extracting metadata, and validating sample alignments using `pandas` and custom helper functions.

---

## Reproducibility 

To ensure computational reproducibility, `.yml` environment files and `config.sh` scripts are provided within the project directories. 

To replicate an environment using Conda:
```bash
conda env create -f scRNAseq/environment.yml
conda activate sc-rna-env

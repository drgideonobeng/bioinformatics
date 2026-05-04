# Single-Cell ATAC-seq Gene Regulatory Network Pipeline

## Overview
This is a complete, optimized 10-step pipeline for analyzing Single-Cell ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) data. Built primarily in R using **Seurat** and **Signac**, this workflow takes raw genomic sequencing fragments and processes them through rigorous quality control, dimensionality reduction, cell-type annotation, and motif analysis, culminating in the construction of unbiased 3D Gene Regulatory Networks (GRNs).

## Key Findings (PBMC Dataset)
* **High-Resolution Clustering:** Successfully partitioned complex Peripheral Blood Mononuclear Cells (PBMCs) into highly distinct biological compartments, notably isolating **CD14+ Monocytes**.
* **Transcription Factor Discovery:** Identified master regulatory motifs—specifically **IRF1**, **Arid3a**, **ZNF354C**, and **ZNF384**—driving the monocyte-specific chromatin landscape.
* **Single-Cell Motif Activity:** Utilized ChromVAR to map transcription factor activity at single-cell resolution, confirming that IRF1 acts as an active master regulator exclusively within the myeloid lineage.
* **Enhancer-Promoter Looping (GRNs):** Conducted unbiased discovery of monocyte target genes (e.g., *RAB31*, *RP11-73M18.2*) and mapped physical 3D enhancer-to-promoter linkages, successfully reverse-engineering the regulatory wiring of the cell.

## Pipeline Architecture

### Part 1: Preprocessing & Clustering
* **Step 01-03: Quality Control** - Filtering fragments based on peak counts, Nucleosome Banding, and Transcriptional Start Site (TSS) enrichment.
* **Step 04: Dimensionality Reduction** - Term Frequency-Inverse Document Frequency (TF-IDF) normalization followed by Singular Value Decomposition (SVD) to perform Latent Semantic Indexing (LSI).
* **Step 05: Graph-Based Clustering** - Non-linear dimension reduction via UMAP to group cells with similar chromatin landscapes.
* **Step 06: Cell Type Annotation** - Inferring Gene Activity from ATAC reads to accurately annotate distinct lineages (e.g., Dendritic cells, B cells, T cells, Monocytes).

### Part 2: Regulatory Discovery & GRNs
* **Step 07: Differential Accessibility & Peak Calling** - Identifying monocyte-specific genomic regions (marker peaks).
* **Step 08: Motif Enrichment** - Scanning differentially accessible peaks against the JASPAR database to identify enriched transcription factor binding sequences.
* **Step 09: ChromVAR Activity Mapping** - Calculating background-corrected, single-cell motif activity scores to visualize master regulators directly on the UMAP.
* **Step 10: Peak-to-Gene Linking** - Dynamically discovering target genes and running cell-by-cell Pearson correlations to construct definitive 3D chromatin loops (Coverage Plots).

## Project Structure
```text
├── scripts/
│   ├── 01_to_06_preprocessing_and_clustering.R
│   ├── 07_differential_accessibility.R
│   ├── 08_motif_enrichment.R
│   ├── 09_chromVAR_activity.R
│   └── 10_peak_to_gene_linking.R
├── results/
│   ├── objects/         # Checkpointed Seurat .rds files
│   ├── tables/          # Discovered target genes and enriched motif CSVs
│   └── plots/           # Output PDFs (UMAPs, QC Violins, Coverage/Looping Plots)
└── README.md

### Usage
Scripts are designed to be run sequentially from the command line. Multicore parallelization (future) is implemented by default to significantly reduce computation time on steps like genome-wide peak-to-gene linkage.

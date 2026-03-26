# Modular Single-Cell RNA-Seq Pipeline

This repository contains a modular, reproducible pipeline for processing single-cell RNA-sequencing (scRNA-seq) data using the `Seurat` package in R. 

It is designed to take raw count matrices and process them through quality control, normalization, dimensionality reduction, clustering, marker discovery, and final cell type annotation.

## Project Structure

The pipeline is heavily organized to keep code, raw data, and generated results separate.

├── config.sh                  # Master configuration file (sets variables & thresholds)
├── data/                      # Raw input data (e.g., Cell Ranger outputs)
├── results/                   
│   ├── objects/               # Saved Seurat objects (.rds) at each major step
│   ├── plots/                 # Generated PDFs (UMAPs, FeaturePlots, etc.)
│   └── tables/                # Exported data (e.g., Marker gene CSVs)
└── scripts/                   # Modular R scripts
    ├── 01_... to 04_...       # Setup, loading, and QC filtering
    ├── 05_normalize_data.R    # Normalization, Variable Features, and Scaling
    ├── 06_run_dim_reduction.R # PCA, Elbow Plot, Neighbors, Clusters, and UMAP
    ├── 07_plot_umap.R         # Visualizing PCA and UMAP clusters
    ├── 08_find_markers.R      # Calculating differential gene expression
    ├── 09_plot_marker_genes.R # Visualizing textbook markers via FeaturePlot
    └── 10_annotate_clusters.R # Assigning biological names to clusters

## How to Use This Pipeline

### 1. Set Your Configuration
Before running any scripts, open `config.sh` to define your input directories, output directories, and biological thresholds (e.g., mitochondrial percentage, PCA dimensions, clustering resolution).

### 2. Load the Environment
Load the variables into your terminal session:
```bash
source config.sh

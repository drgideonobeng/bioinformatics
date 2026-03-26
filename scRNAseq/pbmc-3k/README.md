# Automated scRNA-seq Processing Pipeline: PBMC 3k

An end-to-end, fully automated Single-Cell RNA-Sequencing (scRNA-seq) analysis pipeline built with **Bash** and **R (Seurat)**. This project processes the classic 10x Genomics Peripheral Blood Mononuclear Cells (PBMC) 3k dataset from raw count matrices to biologically annotated cell clusters.



## Biological Overview
Peripheral blood contains a diverse population of immune cells. This pipeline ingests raw transcriptomic data and computationally isolates distinct biological populations, including:
* Naïve and Memory CD4+ T Cells
* CD8+ T Cells
* B Cells
* CD14+ and FCGR3A+ Monocytes
* Natural Killer (NK) Cells
* Dendritic Cells (DCs) and Platelets

##️ Pipeline Architecture
This project utilizes a **Decoupled Architecture**, separating configuration variables from execution scripts to ensure maximum reproducibility and ease of use.

* **`config.sh`**: The single source of truth. Contains all dynamic parameters (thresholds, resolutions, PCA dims).
* **`run_pipeline.sh`**: The master controller. Sources the config and orchestrates the R scripts in sequence.
* **`scripts/`**: Modular, single-purpose R scripts utilizing the `fs` and `tidyverse` ecosystems for robust file routing and data manipulation.

### Directory Structure
```text
pbmc/
├── README.md
├── .gitignore
├── config.sh                   # Environment variables & pipeline parameters
├── run_pipeline.sh             # Master execution script
└── scripts/
    ├── 01_data_download.R      # Fetches 10x Genomics tarball
    ├── 02_create_seurat_obj.R  # Initializes Seurat object & calculates MT%
    ├── 03_qc_visualize.R       # Generates pre-filter violin/scatter plots
    ├── 04_filter_cells.R       # Subsets cells based on config.sh thresholds
    ├── 05_run_dim_reduction.R  # PCA, Neighbor finding, and UMAP
    ├── 06_visualize_clusters.R # Exports UMAP and PCA plots
    ├── 07_find_markers.R       # Differential expression for cluster biomarkers
    └── 08_annotate_cells.R     # Maps biological identities to clusters

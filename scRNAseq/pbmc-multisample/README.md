# Single-Cell RNA-Seq Pipeline: PBMC Multi-Sample Analysis

## Overview
This repository contains a complete, end-to-end Single-Cell RNA-Sequencing (scRNA-seq) analysis pipeline built in R using **Seurat v5**. It was designed to process, integrate, and analyze human Peripheral Blood Mononuclear Cells (PBMCs) across two distinct conditions: **Untreated** vs. **Treated/Stimulated**.

The primary biological finding of this pipeline is the identification of a massive, coordinated **Viral/Interferon Response** (driven by genes like *ISG15*, *IFIT1*, and *IFIT3*) specifically within the B-cell population after exposure to the stimulus.

## Project Structure
```text
pbmc_multi_sample/
├── data/                   # Raw Cell Ranger outputs (Untreated & Treated)
├── results/
│   ├── objects/            # Saved Seurat .rds files (Checkpoints)
│   ├── plots/              # PDF visualizations (UMAPs, Volcano, Pathways)
│   └── tables/             # CSVs of differential genes and GO terms
├── scripts/                # The 11-step analysis pipeline
├── env.yml                 # Conda environment dependencies
└── README.md               # Project documentation

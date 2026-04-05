# Arabidopsis RNA-Seq Quantification: A Salmon Pseudoalignment Pipeline

# Arabidopsis Thaliana RNA-Seq Analysis Pipeline

This project performs a differential expression (DE) analysis on *Arabidopsis thaliana* samples to identify transcriptomic changes under specific experimental conditions (Roots vs. Shoots / Treated vs. Control).

##  Summary of Analysis
- **Reference Genome:** *A. thaliana* (TAIR10)
- **Quantification:** Salmon (mapping-based)
- **Differential Expression:** DESeq2 (R/Bioconductor)
- **Annotation:** Org.At.tair.db / biomaRt

## Project Structure
- `scripts/`: Shell scripts for data download, QC (FastQC), and Salmon quantification.
- `rscripts/`: R scripts for DESeq2 analysis, batch correction, and visualization.
- `metadata.csv`: Sample mapping and experimental design.
- `results/`: Processed output, including:
    - `DE_results_Annotated.csv`: Final list of differentially expressed genes.
    - `plots/`: Visualizations including Volcano plots, PCA, and Heatmaps.

## How to Run
1. Update `config.sh` with your local paths.
2. Run the master pipeline:
   ```bash
   bash run_pipeline.sh
   bash run_r_pipeline.sh

# Bioinformatics and Computational Biology Portfolio

Welcome to my bioinformatics repository. This space serves as a centralized hub for the computational pipelines and workflows I have developed to investigate complex biological questions.

As a trained cell biologist, I approach data with a deep curiosity about how molecular signatures influence cellular fate. I am passionate about bridging the gap between high-throughput sequencing data and meaningful biological insights. This repository showcases end-to-end, reproducible computational pipelines. By combining my extensive understanding of cellular biology with tools such as R, Python, and Bash, I design automated workflows to process high-throughput multi-omic datasets, conduct rigorous statistical analyses, and generate publication-ready visualizations.

## Technical Skills

- **Languages:** R (Bioconductor), Python (Pandas, NumPy), Bash Scripting
- **Single-Cell Analysis:** Seurat (v4/v5), Signac, Trajectory Inference (Slingshot), Harmony Integration
- **Transcriptomics:** DESeq2, Salmon, STAR, clusterProfiler (GO/KEGG Enrichment)
- **Genomics:** GATK Best Practices, Variant Discovery (BWA, HaplotypeCaller)
- **Functional Genomics:** Multi-omic integration and CRISPR-Cas9 screen analysis
- **Infrastructure:** Conda environment management, Git version control, modular pipeline architecture

## Project Highlights

### Single-Cell Transcriptomics and Developmental Modeling

My work focuses on mapping cellular heterogeneity and the logic of lineage commitment.

- **Mesp1 Cohort Analysis:** I developed a trajectory inference pipeline to model murine vascular development from E8.0 to E10.5. By integrating multiple datasets and applying Slingshot pseudotime, I identified the transition from early proliferative progenitors to mature endothelium. These findings were validated by enriching oxygen-response and fluid shear-stress signaling pathways.

- **Perturbation Studies:** I built specialized pipelines for the mesp1CreTbx1KO and MSC-HUVEC datasets to evaluate how specific genetic knockouts or varied environmental conditions alter cell-to-cell communication and local differential expression.

### Single-Cell Epigenomics (scATAC-seq)

I design workflows to uncover the regulatory blueprints behind gene expression. This includes processing raw fragments, performing dimensionality reduction via Latent Semantic Indexing (LSI), and peak calling to identify differential chromatin accessibility across distinct cell states.

### Bulk RNA-Seq and Variant Discovery

- **Cancer Transcriptomics:** Analysis of colon cancer datasets, including driver gene identification and survival probability modeling.

- **GATK Pipelines:** Implementation of the Broad Institute's Best Practices for variant calling, covering everything from raw read alignment to high-confidence SNP and Indel identification in both DNA and RNA sequencing data.

### Functional Genomics and CRISPR Screens

I have automated the quantification and statistical evaluation of high-throughput CRISPR screens. These scripts are designed to identify essential genetic dependencies across different cellular contexts, providing an added functional layer to transcriptomic findings.

## Biological Intuition and Data Rigor

Bioinformatics is more than just running software; it is about ensuring that computational outputs align with biological reality. I prioritize:

- **Biological Validation:** I conduct sanity checks on every trajectory and cluster against known developmental markers and metabolic shifts to ensure that the results are grounded in biology.

- **Scalable Automation:** I utilize configuration scripts and wrapper functions to ensure that pipelines can handle expanding datasets with minimal manual intervention.

- **Reproducible Research:** I document every step through environment files and modular code to ensure that my findings are fully transparent and replicable by other researchers.

## Reproducibility and Usage

To replicate any of the environments used in these projects, please use the provided Conda environment files.

### Bash Instructions

```bash
# Clone the repository
git clone https://github.com/drgideonobeng/bioinformatics.git
cd bioinformatics

# Build and activate the primary single-cell environment
conda env create -f scRNAseq/environment.yml
conda activate sc-rna-env
```

Feel free to reach out if you have any questions or need further information!

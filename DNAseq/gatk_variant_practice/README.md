# GATK Variant Discovery Pipeline (NA12878)

This repository contains a full genomic analysis pipeline using **GATK4** to identify variants in human sample NA12878 (Chromosome 20).

## 🧬 Project Overview
The pipeline processes raw sequence data through alignment, quality recalibration, and variant calling, culminating in the identification of high-confidence SNPs and Indels.

## 🛠 Tech Stack
* **Alignment:** BWA-MEM
* **Processing:** Samtools, Picard
* **Variant Calling:** GATK HaplotypeCaller
* **Filtering:** GATK Hard Filtering (QD, FS, SOR, MQ, etc.)
* **Visualization:** IGV (Integrative Genomics Viewer)

## 📊 Key Results
* **Transition/Transversion (Ts/Tv) Ratio:** 1.94 (indicative of high-quality human data).
* **Significant Finding:** Identified a **Homozygous SNP (T→C)** at `chr20:88,108` within the **DEFB125** (Defensin Beta 125) gene. 
* **Verification:** Confirmed via IGV with an Allele Frequency of 1.0 and a Quality score of 681.06.

## 🚀 How to Run
1. Setup environment: `conda env create -f env/conda_environment.yml`
2. Configure paths in `scripts/config.sh`
3. Execute scripts sequentially from `05_index_reference.sh` to `70_annotate_variants.sh`

## ⚠️ Troubleshooting Note
On some macOS/Rutgers systems, the shell may exit due to `HISTTIMEFORMAT` errors. If the terminal closes unexpectedly, run:
`set +eu`

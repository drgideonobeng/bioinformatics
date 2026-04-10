#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

cat("=== Step 10: Generating Final Research Report ===\n")

# 1. Create a Clean Report Directory
dir_create("final_report/figures")
dir_create("final_report/data")

# 2. Identify and Copy Figures
figures_to_copy <- dir_ls("figures", glob = "*.png")
file_copy(figures_to_copy, "final_report/figures/", overwrite = TRUE)

# 3. Harvest Top 10 Pathways from GSEA Results
if(file_exists("results/gsea_hallmark_results.csv")) {
  gsea <- read_csv("results/gsea_hallmark_results.csv", show_col_types = FALSE)
  
  # clusterProfiler uses 'p.adjust', not 'padj'
  top_10 <- gsea %>%
    arrange(p.adjust) %>% 
    slice_head(n = 10) %>%
    select(ID, NES, p.adjust)
  
  write_csv(top_10, "final_report/data/top_10_pathways.csv")
}

# 4. Generate a Project README Summary
readme_content <- "
# RNA-Seq Analysis: Human Airway Response to Dexamethasone

## Project Overview
This pipeline analyzed the transcriptional response of Human Airway Smooth Muscle (ASM) 
cells to 6.0 uM Dexamethasone treatment.

## Key Findings
- **Differential Expression:** Identified signature genes (e.g., DUSP1, CRISPLD2).
- **Primary Activation:** Adipogenesis and Metabolic remodeling.
- **Primary Suppression:** Inflammatory Response and TNFA Signaling via NFkB.

## Directory Structure
- `/figures`: Volcano, Heatmap, GSEA, and Network plots.
- `/data`: Top pathway statistics.
- `/scripts`: Reproducible R code.
"

write_lines(readme_content, "final_report/README.md")

cat("[OK] Project organized! Check the 'final_report/' folder.\n")

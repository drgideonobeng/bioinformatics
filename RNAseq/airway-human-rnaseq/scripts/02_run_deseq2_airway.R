#!/usr/bin/env Rscript

# --- Step 1: Load Libraries ---
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(fs)
})

cat("=== Step 2: Running DESeq2 Analysis (Airway) ===\n\n")

# --- Step 2: Load Data ---
# Note: Using relative paths assuming we run this from the project root
counts <- read_csv("data/raw/airway_counts.csv", show_col_types = FALSE) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

metadata <- read_csv("metadata/airway_metadata.csv", show_col_types = FALSE) %>%
  column_to_rownames("sample_id")

# --- Step 3: Setup DESeq2 Object ---
# Ensure 'Control' is the reference level
metadata$condition <- factor(metadata$condition, levels = c("Control", "Treated"))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

# --- Step 4: Pre-filtering ---
# We remove genes that have very few reads. This improves statistical power.
# Rule of thumb: Keep genes with at least 10 reads total across all samples.
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("[INFO] Genes remaining after filtering:", nrow(dds), "\n")

# --- Step 5: Run DESeq2 ---
# Because this is real data, the standard 'DESeq()' function will work 
# perfectly without the manual fixes we used for the practice data.
cat("[INFO] Running DESeq() pipeline...\n")
dds <- DESeq(dds)

# --- Step 6: Extract & Shrink Results ---
# We use 'apeglm' shrinkage for high-quality Log Fold Changes.
cat("[INFO] Shrinking results with apeglm...\n")
res <- lfcShrink(dds, coef = "condition_Treated_vs_Control", type = "apeglm")

# --- Step 7: Save Results ---
dir_create("results")

res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  arrange(padj) # Put most significant genes at the top

write_csv(res_df, "results/airway_de_results.csv")
cat("[OK] Analysis complete! Results saved to 'results/airway_de_results.csv'.\n")

# Save the full dds object for future use (like heatmaps)
saveRDS(dds, "results/dds_airway.rds")

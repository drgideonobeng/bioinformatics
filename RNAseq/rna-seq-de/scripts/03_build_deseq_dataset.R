#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
})

cat("=== RNA-seq Step 2C: Build DESeq2 Dataset (simplified) ===\n\n")

# 1) Load data
cat("[INFO] Loading counts + metadata...\n")
counts_df   <- read_csv("data/raw/practice_counts_matrix.csv", show_col_types = FALSE)
metadata_df <- read_csv("metadata/practice_sample_metadata.csv", show_col_types = FALSE)

# 2) Prepare DESeq2 inputs
cat("[INFO] Preparing count matrix and metadata for DESeq2...\n")

count_mat <- counts_df %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round()

coldata <- metadata_df %>%
  column_to_rownames("sample_id") %>%
  mutate(condition = factor(condition))

# 3) Safety check: sample alignment
cat("[INFO] Checking sample alignment...\n")
stopifnot(identical(colnames(count_mat), rownames(coldata)))
cat("[OK] Samples aligned.\n")

# 4) Build DESeq2 dataset
cat("[INFO] Building DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
cat("[OK] DESeqDataSet created.\n\n")

# 5) Quick summary
cat("[INFO] Dataset summary:\n")
print(dds)

cat("\n[INFO] Condition levels:\n")
print(levels(coldata$condition))

cat("\n[OK] Step 2C complete.\n")

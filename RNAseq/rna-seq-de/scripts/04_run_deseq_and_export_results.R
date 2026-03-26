#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(DESeq2)
  library(apeglm)
})

cat("=== RNA-seq Step 2D: Run DESeq2 + Export Results ===\n\n")

# -----------------------------
# 1) Load counts + metadata
# -----------------------------
counts_file <- "data/raw/practice_counts_matrix.csv"
metadata_file <- "metadata/practice_sample_metadata.csv"

cat("[INFO] Loading counts + metadata...\n")
counts_df   <- read_csv(counts_file, show_col_types = FALSE)
metadata_df <- read_csv(metadata_file, show_col_types = FALSE)

# -----------------------------
# 2) Prepare DESeq2 inputs
# -----------------------------
cat("[INFO] Preparing count matrix and metadata for DESeq2...\n")

count_mat <- counts_df %>%
  column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round()

coldata <- metadata_df %>%
  column_to_rownames("sample_id") %>%
  mutate(condition = factor(condition))

cat("[INFO] Checking sample alignment...\n")
stopifnot(identical(colnames(count_mat), rownames(coldata)))
cat("[OK] Samples aligned.\n")

# -----------------------------
# 3) Build DESeq2 dataset
# -----------------------------
cat("[INFO] Building DESeqDataSet...\n")
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)

# Optional: filter very low-count genes (common practice)
dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("[OK] DESeqDataSet ready.\n\n")

# -----------------------------
# 4) Run DESeq2
# -----------------------------
cat("[INFO] Running DESeq() model fitting...\n")
dds <- DESeq(dds)
cat("[OK] DESeq() complete.\n\n")

# -----------------------------
# 5) Extract results (Treatment vs Control)
# -----------------------------
cat("[INFO] Extracting results: Treatment vs Control...\n")
res <- results(dds, contrast = c("condition", "Treatment", "Control"))

res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

cat("[OK] Results extracted.\n\n")

# -----------------------------
# 6) Shrink log2 fold changes (apeglm)
# -----------------------------
cat("[INFO] Shrinking log2 fold changes using apeglm...\n")
res_shrunk <- lfcShrink(dds, coef = "condition_Treatment_vs_Control", type = "apeglm")

res_shrunk_tbl <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id") %>%
  as_tibble()

cat("[OK] Shrinkage complete.\n\n")

# -----------------------------
# 7) Save outputs
# -----------------------------
dir.create("results", showWarnings = FALSE)

out_unshrunk <- "results/deseq2_results_unshrunk.csv"
out_shrunk   <- "results/deseq2_results_shrunk_apeglm.csv"

write_csv(res_tbl, out_unshrunk)
write_csv(res_shrunk_tbl, out_shrunk)

cat("[OK] Saved unshrunk results to: ", out_unshrunk, "\n", sep = "")
cat("[OK] Saved shrunk results to:   ", out_shrunk, "\n\n", sep = "")

# -----------------------------
# 8) Print top results (by adjusted p-value)
# -----------------------------
cat("[INFO] Top genes by adjusted p-value (shrunk results):\n")

top_tbl <- res_shrunk_tbl %>%
  arrange(padj) %>%
  select(gene_id, log2FoldChange, lfcSE, stat, pvalue, padj)

print(top_tbl)

cat("\n[OK] Step 2D complete.\n")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

cat("=== RNA-seq Step 2B: Inspect Counts + Metadata (tidyverse) ===\n\n")

# -----------------------------
# 1) Load files
# -----------------------------
counts_file <- "data/raw/practice_counts_matrix.csv"
metadata_file <- "metadata/practice_sample_metadata.csv"

cat("[INFO] Loading count matrix:", counts_file, "\n")
counts_df <- read_csv(counts_file, show_col_types = FALSE)

cat("[INFO] Loading metadata:    ", metadata_file, "\n")
metadata_df <- read_csv(metadata_file, show_col_types = FALSE)

# -----------------------------
# 2) Basic structure checks
# -----------------------------
cat("\n[INFO] Count matrix dimensions (rows, cols): ", nrow(counts_df), ", ", ncol(counts_df), "\n", sep = "")
cat("[INFO] Metadata dimensions (rows, cols):     ", nrow(metadata_df), ", ", ncol(metadata_df), "\n", sep = "")

cat("\n[INFO] Count matrix columns:\n")
print(colnames(counts_df))

cat("\n[INFO] Metadata columns:\n")
print(colnames(metadata_df))

# -----------------------------
# 3) Extract sample names
# -----------------------------
count_sample_names <- colnames(counts_df) %>%
  setdiff("gene_id")

metadata_sample_names <- metadata_df %>%
  pull(sample_id)

cat("\n[INFO] Sample names from count matrix:\n")
print(count_sample_names)

cat("\n[INFO] Sample names from metadata:\n")
print(metadata_sample_names)

# -----------------------------
# 4) Check exact order match
# -----------------------------
if (identical(count_sample_names, metadata_sample_names)) {
  cat("\n[OK] Sample names are aligned (same values, same order).\n")
} else {
  cat("\n[WARN] Sample names are NOT aligned.\n")

  if (setequal(count_sample_names, metadata_sample_names)) {
    cat("[INFO] Same sample names are present, but order differs.\n")
  } else {
    cat("[INFO] Sample name sets differ.\n")
    
    missing_in_metadata <- setdiff(count_sample_names, metadata_sample_names)
    extra_in_metadata <- setdiff(metadata_sample_names, count_sample_names)

    if (length(missing_in_metadata) > 0) {
      cat("[WARN] Missing in metadata:\n")
      print(missing_in_metadata)
    }
    if (length(extra_in_metadata) > 0) {
      cat("[WARN] Extra in metadata:\n")
      print(extra_in_metadata)
    }
  }
}

# -----------------------------
# 5) Check count columns are numeric
# -----------------------------
cat("\n[INFO] Checking count columns are numeric:\n")

numeric_check_tbl <- tibble(sample_col = count_sample_names) %>%
  mutate(is_numeric = sapply(sample_col, function(x) is.numeric(counts_df[[x]])))

print(numeric_check_tbl)

if (all(numeric_check_tbl$is_numeric)) {
  cat("[OK] All count columns are numeric.\n")
} else {
  cat("[WARN] Some count columns are not numeric.\n")
}

# -----------------------------
# 6) Preview data
# -----------------------------
cat("\n[INFO] Count matrix preview:\n")
print(counts_df)

cat("\n[INFO] Metadata preview:\n")
print(metadata_df)

cat("\n[OK] Inspection complete.\n")

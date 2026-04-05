#!/usr/bin/env Rscript

cat("=== RNA-seq Step 2B: Inspect Counts + Metadata ===\n\n")

# -----------------------------
# 1) Load files
# -----------------------------
counts_file <- "data/raw/practice_counts_matrix.csv"
metadata_file <- "metadata/practice_sample_metadata.csv"

cat("[INFO] Loading count matrix:", counts_file, "\n")
counts_df <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)

cat("[INFO] Loading metadata:    ", metadata_file, "\n")
metadata_df <- read.csv(metadata_file, stringsAsFactors = FALSE)

# -----------------------------
# 2) Basic structure checks
# -----------------------------
cat("\n[INFO] Count matrix dimensions (rows, cols): ", nrow(counts_df), ", ", ncol(counts_df), "\n", sep = "")
cat("[INFO] Metadata dimensions (rows, cols):     ", nrow(metadata_df), ", ", ncol(metadata_df), "\n", sep = "")

cat("\n[INFO] Count matrix column names:\n")
print(colnames(counts_df))

cat("\n[INFO] Metadata columns:\n")
print(colnames(metadata_df))

# -----------------------------
# 3) Extract sample names
# -----------------------------
# Count matrix sample columns = all columns except gene_id
count_sample_names <- colnames(counts_df)[colnames(counts_df) != "gene_id"]

# Metadata sample names from sample_id column
metadata_sample_names <- metadata_df$sample_id

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
  
  # Report differences
  if (setequal(count_sample_names, metadata_sample_names)) {
    cat("[INFO] They contain the same sample names, but the order is different.\n")
  } else {
    cat("[INFO] The sample name sets are different.\n")
    
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
for (sample_col in count_sample_names) {
  is_num <- is.numeric(counts_df[[sample_col]])
  if (is_num) {
    cat("[OK]  ", sample_col, " is numeric\n", sep = "")
  } else {
    cat("[FAIL] ", sample_col, " is NOT numeric\n", sep = "")
  }
}

# -----------------------------
# 6) Preview data
# -----------------------------
cat("\n[INFO] Count matrix preview:\n")
print(head(counts_df))

cat("\n[INFO] Metadata preview:\n")
print(metadata_df)

cat("\n[OK] Inspection complete.\n")

#!/usr/bin/env Rscript

cat("=== RNA-seq Step 2A: Create Practice Counts + Metadata ===\n\n")

# -----------------------------
# 1) Create a small count matrix
#    Rows = genes
#    Columns = samples
# -----------------------------
counts_df <- data.frame(
  gene_id = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH"),
  Ctrl_1 = c(100, 250, 400,  50, 800, 120, 300,  90),
  Ctrl_2 = c(110, 240, 390,  45, 780, 130, 310,  95),
  Ctrl_3 = c( 95, 255, 410,  55, 820, 125, 295, 100),
  Trt_1  = c(300, 245, 390, 180, 810, 115, 305,  88),  # GeneA and GeneD up in treatment
  Trt_2  = c(320, 260, 405, 170, 790, 118, 290,  92),
  Trt_3  = c(310, 250, 395, 190, 805, 122, 298,  91),
  check.names = FALSE
)

# -----------------------------
# 2) Create sample metadata
#    One row per sample
# -----------------------------
metadata_df <- data.frame(
  sample_id = c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Trt_1", "Trt_2", "Trt_3"),
  condition = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
  stringsAsFactors = FALSE
)

# -----------------------------
# 3) Save files
# -----------------------------
counts_out <- "data/raw/practice_counts_matrix.csv"
metadata_out <- "metadata/practice_sample_metadata.csv"

write.csv(counts_df, counts_out, row.names = FALSE)
write.csv(metadata_df, metadata_out, row.names = FALSE)

# -----------------------------
# 4) Print summaries
# -----------------------------
cat("[OK] Count matrix saved to:", counts_out, "\n")
cat("[OK] Metadata saved to:    ", metadata_out, "\n\n")

cat("[INFO] Count matrix preview:\n")
print(counts_df)

cat("\n[INFO] Metadata preview:\n")
print(metadata_df)

cat("\n[OK] Practice data creation complete.\n")

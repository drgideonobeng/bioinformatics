#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(stringr) # For readable string manipulation
})

cat("=== Step 4: Mapping Ensembl IDs to Gene Symbols ===\n")

# 1. Load results
res <- read_csv("results/airway_de_results.csv", show_col_types = FALSE)

# 2. Clean IDs and Map (Tidyverse Style)
res_annotated <- res %>%
  # Remove the version decimal (e.g., .14) from the end of the ID
  mutate(gene_id_clean = str_remove(gene_id, "\\..*$")) %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, keys = gene_id_clean, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
    entrez = mapIds(org.Hs.eg.db, keys = gene_id_clean, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  ) %>%
  # Organize columns and remove the temporary clean ID
  select(gene_id, symbol, entrez, log2FoldChange, padj, everything()) %>%
  select(-gene_id_clean)

# 3. Save
write_csv(res_annotated, "results/airway_de_results_annotated.csv")

cat("[OK] Done! Results saved with human-readable symbols.\n")

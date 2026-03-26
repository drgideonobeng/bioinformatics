#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
})

cat("=== Step 8: Generating Biological Summary Report ===\n")

# 1. Load Data
de_results <- read.csv("results/colon_cancer_annotated.csv")
gsea_results <- read.csv("results/gsea_hallmark_results.csv")

# 2. Identify "Leading Edge" Genes for the #1 Pathway
# We'll pull the genes that drove the E2F Targets enrichment
cat("-> Extracting core drivers for E2F Targets...\n")
e2f_pathway <- gsea_results %>% filter(pathway == "HALLMARK_E2F_TARGETS")
# Split the semicolon-separated string we created in Step 5
driver_genes <- unlist(strsplit(as.character(e2f_pathway$leadingEdge), ";"))

# 3. Filter DE results to show ONLY these driver genes
# We want to see how much these specific genes changed
summary_table <- de_results %>%
  filter(symbol %in% driver_genes) %>%
  select(symbol, log2FoldChange, padj, baseMean) %>%
  arrange(desc(log2FoldChange)) %>%
  head(20)

# 4. Save a clean summary table
write.csv(summary_table, "results/top_driver_genes_summary.csv", row.names = FALSE)

# 5. Print a "Executive Summary" to the terminal
cat("\n--- EXECUTIVE BIOLOGICAL SUMMARY ---\n")
cat("Condition: TCGA Colorectal Adenocarcinoma (Tumor vs Normal)\n")
cat(paste("Top Pathway: HALLMARK_E2F_TARGETS (NES:", round(e2f_pathway$NES, 2), ")\n"))
cat("\nTop 5 Driver Genes in this Pathway:\n")
print(head(summary_table, 5))

cat("\n[OK] Summary generated in results/top_driver_genes_summary.csv\n")

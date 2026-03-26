#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(fs)
})

cat("=== Step 8: Generating Gene-Concept Network ===\n")

# 1. Prepare data
res <- read_csv("results/airway_de_results_annotated.csv", show_col_types = FALSE)
gene_list_df <- res %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(log2FoldChange))

gene_list <- gene_list_df$log2FoldChange
names(gene_list) <- gene_list_df$symbol

# 2. Run GSEA 
h_df <- msigdbr(species = "Homo sapiens", collection = "H") %>% 
  select(gs_name, gene_symbol)

gsea_res <- GSEA(gene_list, TERM2GENE = h_df, pvalueCutoff = 0.05, verbose = FALSE)

# 3. Create the Cnetplot (Simplified Arguments)
cat("[INFO] Building network of top pathways...\n")

# Use a basic call first to ensure compatibility
# foldChange is used to color the gene nodes (Red = Up, Blue = Down)
plot <- cnetplot(gsea_res, 
                 showCategory = 10, 
                 foldChange = gene_list) +
  ggtitle("Gene-Pathway Connections (Top 3 Hallmark Sets)")

# 4. Save with high resolution
dir_create("figures")
ggsave("figures/airway_cnetplot.png", plot, width = 12, height = 10, bg = "white")

cat("[OK] Network plot saved to figures/airway_cnetplot.png\n")

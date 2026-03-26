#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(fs)
})

cat("=== Step 9: Generating Enrichment Map (emapplot) ===\n")

# 1. Load data
res <- read_csv("results/airway_de_results_annotated.csv", show_col_types = FALSE)
gene_list_df <- res %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(log2FoldChange))

gene_list <- gene_list_df$log2FoldChange
names(gene_list) <- gene_list_df$symbol

# 2. Run GSEA with a relaxed p-value to ensure we get a "Network"
h_df <- msigdbr(species = "Homo sapiens", collection = "H") %>% 
  select(gs_name, gene_symbol)

# Relaxing pvalueCutoff to 0.2 to see the "neighborhood"
gsea_res <- GSEA(gene_list, TERM2GENE = h_df, pvalueCutoff = 0.2, verbose = FALSE)

# 3. Calculate Similarity
cat("[INFO] Calculating similarity between pathways...\n")
gsea_sim <- pairwise_termsim(gsea_res)

# 4. Create the Emapplot
cat("[INFO] Rendering map...\n")
# If the plot fails, we catch the error and print a helpful message
tryCatch({
  plot <- emapplot(gsea_sim, 
                   showCategory = 30, 
                   layout = "nicely") + 
    ggtitle("Enrichment Map: Pathway Connectivity") +
    theme_light()

  # 5. Save
  dir_create("figures")
  ggsave("figures/airway_emapplot.png", plot, width = 14, height = 12, units = "in", bg = "white")
  cat("[OK] Enrichment map saved to figures/airway_emapplot.png\n")
}, error = function(e) {
  cat("[ERROR] Could not generate map. Likely too few pathways connected.\n")
  print(e)
})

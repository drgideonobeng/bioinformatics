#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fgsea)
  library(ggplot2)
  library(tidyverse)
  library(msigdbr)
})

cat("=== Step 7: Plotting Enrichment for E2F Targets ===\n")

# 1. Reload the Rank list (same logic as Step 5)
res_df <- read.csv("results/colon_cancer_annotated.csv") %>% 
  filter(!is.na(symbol) & !is.na(stat)) %>%
  distinct(symbol, .keep_all = TRUE)
ranks <- setNames(res_df$stat, res_df$symbol)

# 2. Get the Pathway genes
h_df <- msigdbr(species = "human", category = "H")
h_list <- split(x = h_df$gene_symbol, f = h_df$gs_name)

# 3. Plot the specific pathway
plot <- plotEnrichment(h_list[["HALLMARK_E2F_TARGETS"]], ranks) +
  labs(title = "Enrichment Plot: HALLMARK_E2F_TARGETS",
       subtitle = "Strong upregulation in TCGA-COAD Tumors")

# 4. Save
ggsave("plots/gsea_enrichment_E2F.png", plot, width = 7, height = 5, bg = "white")

cat("[OK] Enrichment plot saved to plots/gsea_enrichment_E2F.png\n")

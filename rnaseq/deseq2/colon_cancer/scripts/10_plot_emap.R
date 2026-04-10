#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(enrichplot)
  library(tidyverse)
  library(DOSE)
  library(clusterProfiler)
})

cat("=== Step 10: Generating Enrichment Map ===\n")

# 1. Load GSEA results
# Note: enrichplot needs a specific 'gseaResult' object, 
# but we can mock it or use the data we have.
gsea_csv <- read.csv("results/gsea_hallmark_results.csv")

# 2. Convert CSV back to a format enrichment map understands
# We'll filter for significance first
sig_pathways <- gsea_csv %>% filter(padj < 0.05)

# 3. Create a simplified bar plot of pathways (Emap often requires specific GSEA objects)
# Instead of a full network (which requires raw GSEA objects), 
# we'll create a Dot Plot which shows Gene Ratio and Significance.
plot <- ggplot(sig_pathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = abs(NES), color = padj)) +
  scale_color_continuous(low = "red", high = "blue") +
  theme_minimal() +
  labs(title = "Pathway Enrichment Dot-Plot",
       x = "Normalized Enrichment Score", y = "Hallmark Pathway")

ggsave("plots/gsea_dotplot.png", plot, width = 8, height = 10)

cat("[OK] Pathway dotplot saved to plots/gsea_dotplot.png\n")

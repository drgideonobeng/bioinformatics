#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

cat("=== Step 6: Plotting GSEA Summary ===\n")

# 1. Load results
input_file <- "results/gsea_hallmark_results.csv"
if (!file.exists(input_file)) stop("Run Step 5 first!")

gsea_res <- read.csv(input_file) %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

# 2. Create the Bar Plot
plot <- ggplot(gsea_res, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "firebrick3", "FALSE" = "dodgerblue3"), 
                    labels = c("Downregulated", "Upregulated")) +
  theme_minimal() +
  labs(title = "GSEA Hallmark Pathways: TCGA-COAD",
       subtitle = "Tumor vs Normal (padj < 0.05)",
       x = "Pathway", y = "Normalized Enrichment Score (NES)",
       fill = "Direction")

# 3. Save
if(!dir.exists("plots")) dir.create("plots")
ggsave("plots/gsea_summary_bar.png", plot, width = 10, height = 8, bg = "white")

cat("[OK] Summary plot saved to plots/gsea_summary_bar.png\n")

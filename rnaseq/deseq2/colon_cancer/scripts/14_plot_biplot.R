#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(tidyverse)
  library(ggrepel)
})

cat("=== Step 14: Generating PCA Biplot (Samples + Gene Loadings) ===\n")

# 1. Load Data
dds <- readRDS("results/dds_output.rds")
vsd <- vst(dds, blind = FALSE)
pc1_drivers <- read.csv("results/pc1_driver_genes.csv")

# 2. Extract PCA coordinates for samples
rv <- rowVars(assay(vsd))
select_genes <- order(rv, decreasing = TRUE)[seq_len(500)]
pca_res <- prcomp(t(assay(vsd)[select_genes, ]))

# Sample Scores (Dots)
sample_df <- as.data.frame(pca_res$x) %>%
  rownames_to_column("sample_id") %>%
  mutate(condition = colData(dds)$condition)

# Gene Loadings (Arrows)
# We scale the loadings so they fit on the same axes as the samples
load_df <- as.data.frame(pca_res$rotation) %>%
  rownames_to_column("ensembl_id") %>%
  inner_join(pc1_drivers %>% select(ensembl_id, symbol), by = "ensembl_id")

# Multiplier to make arrows visible on the sample scale
mult <- min(
  (max(sample_df$PC1) - min(sample_df$PC1)) / (max(load_df$PC1) - min(load_df$PC1)),
  (max(sample_df$PC2) - min(sample_df$PC2)) / (max(load_df$PC2) - min(load_df$PC2))
) * 0.7

load_df <- load_df %>%
  mutate(x_end = PC1 * mult, y_end = PC2 * mult)

# 3. Plotting
plot <- ggplot() +
  # Draw Samples
  geom_point(data = sample_df, aes(x = PC1, y = PC2, color = condition), alpha = 0.4, size = 2) +
  # Draw Arrows for top 10 most influential genes
  geom_segment(data = load_df %>% head(10), 
               aes(x = 0, y = 0, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size = 0.8) +
  # Label Arrows
  geom_text_repel(data = load_df %>% head(10), 
                  aes(x = x_end, y = y_end, label = symbol),
                  color = "black", fontface = "bold", size = 4) +
  scale_color_manual(values = c("Tumor" = "firebrick3", "Normal" = "dodgerblue3")) +
  theme_minimal() +
  labs(title = "PCA Biplot: Gene-Sample Relationships",
       subtitle = "Arrows indicate genes driving the 24% Variance in PC1",
       x = "PC1 (24%)", y = "PC2 (12%)")

# 4. Save
ggsave("plots/pca_biplot.png", plot, width = 10, height = 8, bg = "white")

cat("[OK] Biplot saved to plots/pca_biplot.png\n")

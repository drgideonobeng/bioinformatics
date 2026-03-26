#!/usr/bin/env Rscript

# scripts/04_volcano_plot.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

cat("=== Step 4: Generating Volcano Plot ===\n")

# 1. Load Annotated Data
res_df <- read.csv("results/colon_cancer_annotated.csv", row.names = 1) %>%
  drop_na(padj)

# 2. Define Significance (Thresholds: Fold Change > 4, Adjusted P < 1e-10)
res_df <- res_df %>%
  mutate(change = case_when(
    padj < 1e-10 & log2FoldChange > 2 ~ "UP",
    padj < 1e-10 & log2FoldChange < -2 ~ "DOWN",
    TRUE ~ "NS"
  ))

# 3. Select top 15 most significant genes to label
top_hits <- res_df %>%
  arrange(padj) %>%
  head(15)

# 4. Plot
plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = change)) +
  geom_point(alpha = 0.4, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "dodgerblue3", "UP" = "firebrick3", "NS" = "grey70")) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(1e-10), linetype = "dashed", alpha = 0.5) +
  # Use the 'symbol' column for labels
  geom_text_repel(data = top_hits, aes(label = symbol), 
                  color = "black", fontface = "bold", size = 4) +
  theme_classic() +
  labs(title = "TCGA-COAD: Gene Expression Profile",
       subtitle = "Tumor vs Normal Samples",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none")

# 5. Save
if(!dir.exists("plots")) dir.create("plots")
ggsave("plots/volcano_plot_final.png", plot, width = 8, height = 6, dpi = 300)

cat("[OK] Volcano plot generated in plots/volcano_plot_final.png\n")

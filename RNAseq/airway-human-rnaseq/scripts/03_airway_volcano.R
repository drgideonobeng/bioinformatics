#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(ggrepel)
})

cat("=== Step 3: Generating Volcano Plot ===\n")

# 1. Load results
res <- read_csv("results/airway_de_results.csv", show_col_types = FALSE)

# 2. Add significance categories for coloring
res <- res %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# 3. Get the Top 10 genes to label them
# We filter for genes that actually have a name/ID and are significant
top_genes <- res %>%
  filter(sig != "Not Significant") %>%
  slice_min(order_by = padj, n = 10)

# 4. Plot
plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Down-regulated" = "steelblue", 
                                "Not Significant" = "grey", 
                                "Up-regulated" = "firebrick")) +
  # Add the labels for top genes
  geom_text_repel(data = top_genes, aes(label = gene_id), 
                  size = 3, color = "black", max.overlaps = 15) +
  theme_minimal() +
  labs(title = "Differential Expression: Human Airway Cells + Dexamethasone",
       subtitle = "Labeled: Top 10 most significant genes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Status") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)

# 5. Create folder if it does not exit
dir_create("results")
dir_create("figures")

# 6. Save the plot
ggsave("figures/airway_volcano_plot.png", plot, width = 8, height = 6, dpi = 300)

cat("[OK] Volcano plot saved to figures/airway_volcano_plot.png\n")

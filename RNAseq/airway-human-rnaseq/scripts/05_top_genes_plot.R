#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
})

cat("=== Step 5: Plotting Top 20 Up-regulated Genes ===\n")

# 1. Load annotated results
res <- read_csv("results/airway_de_results_annotated.csv", show_col_types = FALSE)

# 2. Filter for significant up-regulated genes and take the top 20 by Fold Change
top_20 <- res %>%
  filter(padj < 0.05, !is.na(symbol)) %>%
  slice_max(order_by = log2FoldChange, n = 20)

# 3. Create the bar plot
plot <- ggplot(top_20, aes(x = reorder(symbol, log2FoldChange), y = log2FoldChange)) +
  geom_col(fill = "firebrick") +
  coord_flip() +  # Flip to make gene names easier to read
  theme_minimal() +
  labs(
    title = "Top 20 Genes Up-regulated by Dexamethasone",
    subtitle = "Human Airway Smooth Muscle Cells",
    x = "Gene Symbol",
    y = "Log2 Fold Change (Treated vs Control)"
  )

# 4. Save using the fs package
dir_create("figures")
ggsave("figures/top_20_genes.png", plot, width = 8, height = 7, dpi = 300)

cat("[OK] Bar chart saved to figures/top_20_genes.png\n")

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2) # The standard for high-quality figures
})

cat("=== RNA-seq Step 4: Generate Volcano Plot ===\n\n")

# 1. Load the results we just created
res_df <- read_csv("results/de_results.csv", show_col_types = FALSE)

# 2. Add a column to label genes as 'Significant' or 'Not'
# We'll use the standard pharma threshold: padj < 0.05 and |LFC| > 1
res_df <- res_df %>%
  mutate(status = case_when(
    padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
    padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
    TRUE ~ "Not Significant"
  ))

# 3. Create the plot
cat("[INFO] Creating plot...\n")
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Treatment vs Control",
       subtitle = "Significant genes highlighted in red/blue",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

# 4. Save the figure
ggsave("figures/volcano_plot.png", plot = volcano, width = 8, height = 6)
cat("[OK] Figure saved to 'figures/volcano_plot.png'.\n")

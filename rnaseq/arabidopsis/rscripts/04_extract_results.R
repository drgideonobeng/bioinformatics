#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 4: Differential Expression Analysis (CLE9 vs Water) ---")

# 1. Load your fully analyzed dds object
dds <- readRDS(path("results", "dds_analyzed.rds"))
res_dir <- path("results", "DE_results")
plot_dir <- path("results", "plots")
dir_create(res_dir)
dir_create(plot_dir)

# 2. Extract Results
# We explicitly define the contrast: c("Factor", "Numerator", "Denominator")
# Since we want to see what happens IN the mutant compared TO the wild type:
message("Extracting results for CLE9 vs Water...")
res_unshrunken <- results(dds, contrast = c("Treatment", "CLE9", "Water"), alpha = 0.05)

# 3. Apply Log2 Fold Change Shrinkage
# This is crucial for ranking genes and visualizing them accurately.
message("Applying LFC shrinkage (apeglm)...")
coef_name <- "Treatment_CLE9_vs_Water"
res_shrunken <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res_unshrunken)

# 4. Format and Annotate the Results Table
# Convert to a dataframe, move the rownames (Gene IDs) to an actual column, and sort by P-value
res_df <- as.data.frame(res_shrunken) %>%
  rownames_to_column(var = "GeneID") %>%
  arrange(padj) 

# 5. Filter for Significant DEGs 
# Standard thresholds: Adjusted P-value < 0.05 AND an absolute Log2 Fold Change > 1 (i.e., doubled or halved)
sig_degs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

message(glue("\nSUCCESS: Found {nrow(sig_degs)} significant differentially expressed genes!"))

message("\n--- Extracting Top 20 Genes ---")

# Sort by the absolute value of LFC (to get both massive up- and down-regulated genes)
top_20_genes <- sig_degs %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = 20) %>%
  select(GeneID, log2FoldChange, padj) # Keep only the most important columns

message("The Top 20 Most Drastically Changed Genes:")
print(top_20_genes)

# Save to a quick-reference CSV
write_csv(top_20_genes, path(res_dir, "CLE9_top_20_genes.csv"))

# 6. Save Results to CSV
full_csv <- path(res_dir, "CLE9_vs_Water_all_genes.csv")
sig_csv <- path(res_dir, "CLE9_vs_Water_sig_DEGs.csv")

write_csv(res_df, full_csv)
write_csv(sig_degs, sig_csv)
message(glue("\n=> Saved full results to: {full_csv}"))
message(glue("=> Saved significant DEGs to: {sig_csv}"))

# 7. Visualization: MA Plot
message("\nGenerating MA Plot...")
png(path(plot_dir, "ma_plot_CLE9_vs_Water.png"), width = 800, height = 600, res = 120)
plotMA(res_shrunken, main = "MA Plot: CLE9 vs Water", ylim = c(-5, 5))
dev.off()

# 8. Visualization: Volcano Plot
message("Generating Volcano Plot...")
# Drop NAs for plotting
plot_data <- res_df %>% drop_na(padj, log2FoldChange)

volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot: CLE9 vs Water",
       subtitle = paste("Significant DEGs:", nrow(sig_degs)),
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Significant (padj < 0.05 & |LFC| > 1)") +
  theme(legend.position = "bottom") +
  # Add threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5)

ggsave(filename = path(plot_dir, "volcano_CLE9_vs_Water.png"), plot = volcano_plot, width = 7, height = 6)
message("=> Saved MA and Volcano plots to results/plots/")







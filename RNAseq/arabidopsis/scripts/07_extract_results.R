#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
  library(apeglm)
})

message("--- Step 4: Results & Statistical Visualizations ---")

dds <- readRDS(path("results", "dds_analyzed.rds"))
out_dir <- path("results")
plot_dir <- path("results", "plots")

message("Extracting Treated vs Control contrast...")
res <- results(dds, contrast = c("condition", "Treated", "Control"), alpha = 0.05)

message("\nShowing summary of significant Differentially Expressed Genes")
summary(res)

# Plot 1: MA Plot (Unshrunken)
message("Generating MA Plot...")
png(path(plot_dir, "de_ma_plot.png"), width = 800, height = 600, res = 120)
plotMA(res, main = "MA Plot: Treated vs Control", ylim = c(-8, 8))
dev.off()

# Plot 1a: MA Plot with LFC shrinkage
message("Generating MA Plot with LFC shrinkage to improve the estimates...")
# NOTE: Using 'coef' instead of 'contrast' for the modern apeglm method
coef_name <- "condition_Treated_vs_Control"
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

png(path(plot_dir, "de_ma_lfc_plot.png"), width = 800, height = 600, res = 120)
plotMA(res_shrunk, main = "MA Plot with LFC Shrinkage: Treated vs Control", ylim = c(-5, 5))
dev.off()

# Format and Save CSV (Using the newly shrunk data!)
message("Formatting results table...")
res_tidy <- as.data.frame(res_shrunk) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  drop_na(padj, log2FoldChange) %>%
  arrange(padj)

write_csv(res_tidy, path(out_dir, "DE_results_Treated_vs_Control.csv"))

# Plot 2: Volcano Plot
message("Generating Volcano Plot...")
p_cutoff <- 0.05
fc_cutoff <- 1.0

res_plot <- res_tidy %>%
  mutate(
    Significance = case_when(
      log2FoldChange > fc_cutoff & padj < p_cutoff ~ "Upregulated",
      log2FoldChange < -fc_cutoff & padj < p_cutoff ~ "Downregulated",
      TRUE ~ "Not Sig"
    ),
    Significance = factor(Significance, levels = c("Downregulated", "Not Sig", "Upregulated"))
  )

volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "blue", "Not Sig" = "grey", "Upregulated" = "red")) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(p_cutoff), col = "black", linetype = "dashed") +
  labs(title = "Volcano Plot (Shrunken LFC)", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(filename = path(plot_dir, "de_volcano_plot.png"), plot = volcano, width = 8, height = 6, dpi = 300)

message("======================================================")
message("Pipeline Complete! Differential Expression is finished.")
message("======================================================")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 4: Differential Expression Analysis (Dynamic Contrasts) ---")

out_dir <- Sys.getenv("OUT_DIR")
dds <- readRDS(path(out_dir, "dds_analyzed.rds"))
res_dir <- path(out_dir, "DE_results")
plot_dir <- path(out_dir, "plots")
dir_create(res_dir)

# --- DYNAMIC CONTRAST GENERATION ---
# Because we used relevel() in Step 1, the Control is always index 1, and Treatment is index 2
control_level <- levels(dds$Treatment)[1]
trt_level <- levels(dds$Treatment)[2]

message(glue("Extracting results for {trt_level} vs {control_level}..."))
res_unshrunken <- results(dds, contrast = c("Treatment", trt_level, control_level), alpha = 0.05)

message("Applying LFC shrinkage (apeglm)...")
coef_name <- paste0("Treatment_", trt_level, "_vs_", control_level)
res_shrunken <- lfcShrink(dds, coef = coef_name, type = "apeglm", res = res_unshrunken)

# Format and Annotate
res_df <- as.data.frame(res_shrunken) %>%
  rownames_to_column(var = "GeneID") %>%
  arrange(padj) 

sig_degs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

message(glue("\nSUCCESS: Found {nrow(sig_degs)} significant differentially expressed genes!"))

# Dynamic File Naming based on variables!
prefix <- paste0(trt_level, "_vs_", control_level)
full_csv <- path(res_dir, paste0(prefix, "_all_genes.csv"))
sig_csv <- path(res_dir, paste0(prefix, "_sig_DEGs.csv"))

write_csv(res_df, full_csv)
write_csv(sig_degs, sig_csv)
message(glue("=> Saved outputs to: {res_dir}/"))

# Visualizations
message("\nGenerating MA and Volcano Plots...")
png(path(plot_dir, paste0("ma_plot_", prefix, ".png")), width = 800, height = 600, res = 120)
plotMA(res_shrunken, main = glue("MA Plot: {trt_level} vs {control_level}"), ylim = c(-5, 5))
dev.off()

plot_data <- res_df %>% drop_na(padj, log2FoldChange)
volcano_plot <- ggplot(plot_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05 & abs(log2FoldChange) > 1), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = glue("Volcano Plot: {trt_level} vs {control_level}"),
       subtitle = paste("Significant DEGs:", nrow(sig_degs)),
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Significant") +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5)

ggsave(filename = path(plot_dir, paste0("volcano_", prefix, ".png")), plot = volcano_plot, width = 7, height = 6)






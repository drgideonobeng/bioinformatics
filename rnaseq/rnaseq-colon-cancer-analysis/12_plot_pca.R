#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(tidyverse)
})

cat("=== Step 12: Principal Component Analysis (PCA) ===\n")

# 1. Load the DESeq Object from Step 2
# Note: It is best to use the RDS object saved during the DESeq2 run
if (!file.exists("results/dds_output.rds")) {
  stop("Missing dds_output.rds! Please ensure Step 2 saves the dds object.")
}
dds <- readRDS("results/dds_output.rds")

# 2. Variance Stabilizing Transformation (VST)
# VST is much faster than rlog for 500+ samples
cat("-> Applying Variance Stabilizing Transformation (VST)...\n")
vsd <- vst(dds, blind = FALSE)

# 3. Calculate PCA Data
cat("-> Calculating PCA components...\n")
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# 4. Create the Plot
# We use alpha (transparency) because 504 tumor points will overlap
plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_color_manual(values = c("Tumor" = "firebrick3", "Normal" = "dodgerblue3")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "PCA: TCGA-COAD (n=546)",
       subtitle = "Global Transcriptomic Separation: Tumor vs Normal")

# 5. Save
if(!dir.exists("plots")) dir.create("plots")
ggsave("plots/pca_plot.png", plot, width = 8, height = 6, bg = "white")

cat("[OK] PCA plot saved to plots/pca_plot.png\n")

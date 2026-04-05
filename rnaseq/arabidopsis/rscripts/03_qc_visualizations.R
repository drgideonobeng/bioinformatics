#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(fs)
  library(glue)
})

message("--- Step 3: Quality Control Visualizations ---")

dds <- readRDS(path("results", "dds_analyzed.rds"))
vsd <- readRDS(path("results", "vsd_transformed.rds"))
plot_dir <- path("results", "plots")
dir_create(plot_dir)

# Plot 1: Dispersion
message("Generating Dispersion Plot...")
png(path(plot_dir, "qc_dispersion.png"), width = 800, height = 600, res = 120)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Plot 2: PCA
message("Generating PCA Plot...")
# We now group strictly by Treatment instead of Genotype/Batch
pca_data <- plotPCA(vsd, intgroup = "Treatment", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 5, alpha = 0.9) +
  theme_minimal() +
  labs(title = "PCA: Samples by Treatment",
       x = glue("PC1: {percentVar[1]}% variance"),
       y = glue("PC2: {percentVar[2]}% variance"))

ggsave(filename = path(plot_dir, "qc_pca_plot.png"), plot = pca_plot, width = 7, height = 5)

# Plot 3: Sample Distance Heatmap
message("Generating Sample Distance Heatmap...")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Clean row names to show only Treatment and Sample ID (e.g., "CLE9-SRR37530274")
rownames(sampleDistMatrix) <- paste(vsd$Treatment, rownames(sampleDistMatrix), sep="-")
colnames(sampleDistMatrix) <- NULL

png(path(plot_dir, "qc_sample_heatmap.png"), width = 800, height = 600, res = 120)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         main = "Sample Distances")
dev.off()

message("=> All QC plots saved to results/plots/")

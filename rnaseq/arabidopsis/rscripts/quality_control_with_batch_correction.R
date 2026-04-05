#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(fs)
  library(glue)
  library(limma) # Required for removeBatchEffect
})

message("--- Step 3: Quality Control Visualizations (Batch Corrected) ---")

dds <- readRDS(path("results", "dds_analyzed.rds"))
vsd <- readRDS(path("results", "vsd_transformed.rds"))
plot_dir <- path("results", "plots")
dir_create(plot_dir)

# --- THE FIX: Visually remove the batch effect from the VSD object ---
message("\nApplying Limma Batch Effect Correction to VST data for visualization...")
# Note: We only do this for the VSD object used in plotting. DESeq2 handles the 
# batch effect mathematically internally for the actual differential expression testing.
assay(vsd) <- removeBatchEffect(assay(vsd), batch = vsd$Batch)

# Plot 1: Dispersion
message("Generating Dispersion Plot...")
png(path(plot_dir, "qc_dispersion.png"), width = 800, height = 600, res = 120)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Plot 2: PCA
message("Generating PCA Plot...")
pca_data <- plotPCA(vsd, intgroup = c("Genotype", "Batch"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Genotype, shape = Batch)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA: Samples by Genotype (Batch Effect Removed)",
       x = glue("PC1: {percentVar[1]}% variance"),
       y = glue("PC2: {percentVar[2]}% variance"))

ggsave(filename = path(plot_dir, "qc_pca_plot.png"), plot = pca_plot, width = 7, height = 5)

# Plot 3: Sample Distance Heatmap
message("Generating Sample Distance Heatmap...")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Clean row names to show Genotype, Batch, and Sample ID
rownames(sampleDistMatrix) <- paste(vsd$Genotype, paste0("Batch", vsd$Batch), rownames(sampleDistMatrix), sep="-")
colnames(sampleDistMatrix) <- NULL

png(path(plot_dir, "qc_sample_heatmap.png"), width = 800, height = 600, res = 120)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, 
         main = "Sample Distances (Batch Corrected)")
dev.off()

message("=> All batch-corrected QC plots saved to results/plots/")

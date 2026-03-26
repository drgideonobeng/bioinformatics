# scripts/06_visualize_clusters.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "03_pbmc_clustered.rds")
message("Loading clustered data from: ", load_path)
pbmc <- readRDS(load_path)

# 3. Generate Dimensionality Reduction Plots
message("Generating PCA and UMAP plots...")

# DimPlot is Seurat's built-in way to plot dimensionality reductions
p1 <- DimPlot(pbmc, reduction = "pca") + 
  ggtitle("PCA of PBMC Clusters") +
  theme_bw()

p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + 
  ggtitle("UMAP of PBMC Clusters") +
  theme_bw()

# 4. Save to PDF (One per page)
dir_create(plot_dir)
pdf_path <- path(plot_dir, "02_clustering_plots.pdf")

pdf(pdf_path, width = 8, height = 6)
print(p1)
print(p2)
dev.off()

message("Success! Clustering plots saved to: ", pdf_path)





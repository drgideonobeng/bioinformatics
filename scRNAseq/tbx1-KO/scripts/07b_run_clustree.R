#!/usr/bin/env Rscript

# Install clustree if you don't have it yet
if (!requireNamespace("clustree", quietly = TRUE)) {
  install.packages("clustree", repos = "http://cran.us.r-project.org")
}

library(Seurat)
library(clustree)
library(ggplot2)
library(fs)

# 1. Setup Paths
obj_dir  <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load the object right before or after your initial clustering
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message("Loading Seurat object...")
seurat_obj <- readRDS(load_path)

# 3. Calculate a gradient of resolutions (0.2 up to 1.2)
message("Testing multiple resolutions...")
# Note: Seurat saves each of these as a metadata column named e.g., 'RNA_snn_res.0.2'
seurat_obj <- FindClusters(seurat_obj, resolution = seq(0.2, 1.2, by = 0.2))

# 4. Generate the Clustree Visualization
message("Drawing the cluster tree...")
# The prefix tells clustree which metadata columns to look for
tree_plot <- clustree(seurat_obj, prefix = "RNA_snn_res.") +
  ggtitle("Clustering Resolution Stability Tree") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 5. Save the Plot
dir_create(plot_dir)
pdf_path <- path(plot_dir, "07b_clustree_resolutions.pdf")
pdf(pdf_path, width = 10, height = 12) 
print(tree_plot)
invisible(dev.off())

message("Clustree visualization saved to: ", pdf_path)

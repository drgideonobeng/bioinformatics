# scripts/07_visualize_clusters.R
library(Seurat)
library(fs)
library(ggplot2) # For tweaking the plot title

# 1. Pull dynamic variables from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "04_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Generate the UMAP plot
message("Generating UMAP plot...")
# DimPlot automatically colors by cluster. We add labels so we can see the cluster numbers!
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP: 1,024 Mouse Brain Cells")

# 4. Save the plot to a PDF
dir_create(plot_dir) # Just in case the folder doesn't exist
save_path <- path(plot_dir, "03_umap_clusters.pdf")

pdf(save_path, width = 8, height = 6)
print(umap_plot)
# R quirk: null device just means the PDF saved successfully!
dev.off() 

message("Step 07 Complete. UMAP plot saved to: ", save_path)


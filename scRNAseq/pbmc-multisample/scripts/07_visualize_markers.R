# scripts/07_visualize_markers.R
library(Seurat)
library(ggplot2)
library(fs)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 07...")
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

# ========== 2. LOAD UMAP OBJECT ==========
message("Loading clustered UMAP dataset...")
umap_obj <- readRDS(path(obj_dir, "05_umap.rds"))

message("Joining Seurat v5 data layers for visualization...")
umap_obj <- JoinLayers(umap_obj)

# ========== 3. DEFINE MARKERS TO PLOT ==========
message("Defining key biomarker genes...")
# REMOVED the duplicate CD79A
genes_to_plot <- c("GIMAP7", "CCL8", "CD79A", "GNLY", "MS4A4A", "CLIC3", "VMO1", "SELL")

# ========== 4. GENERATE SPLIT UMAP (The 5 vs 6 Answer!) ==========
message("Generating UMAP split by condition...")
pdf(path(plot_dir, "07a_Condition_UMAP.pdf"), width = 12, height = 6)
# orig.ident contains the folder names (untreated vs treated)
print(DimPlot(umap_obj, reduction = "umap", split.by = "orig.ident", label = TRUE))
dev.off()

# ========== 5. GENERATE FEATURE PLOTS ==========
message("Generating FeaturePlots (coloring cells by gene expression)...")
pdf(path(plot_dir, "07b_FeaturePlots.pdf"), width = 12, height = 8)
print(FeaturePlot(umap_obj, features = genes_to_plot, ncol = 3, pt.size = 0.5))
dev.off()

# ========== 6. GENERATE DOT PLOT ==========
message("Generating DotPlot (gene expression across all clusters)...")
pdf(path(plot_dir, "07c_DotPlot.pdf"), width = 10, height = 6)
print(DotPlot(umap_obj, features = genes_to_plot) + RotatedAxis())
dev.off()

message("Step 07 Complete! FeaturePlots, DotPlots, and Split UMAPs saved.")

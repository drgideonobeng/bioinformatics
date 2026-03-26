# scripts/05_cluster_and_umap.R
library(Seurat)
library(ggplot2)
library(fs)
library(glue)
library(patchwork)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 05...")
obj_dir     <- Sys.getenv("OBJ_DIR")
plot_dir    <- Sys.getenv("PLOT_DIR")

pca_dims    <- as.numeric(Sys.getenv("PCA_DIMS"))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES"))
umap_title  <- Sys.getenv("UMAP_TITLE")

# ========== 2. LOAD PCA OBJECT ==========
message("Loading PCA-reduced dataset...")
pca_obj <- readRDS(path(obj_dir, "04_pca.rds"))

# ========== 3. CLUSTER THE CELLS ==========
message(glue("Building Shared Nearest Neighbor (SNN) graph using {pca_dims} dimensions..."))
# This connects cells that are biologically similar
cluster_obj <- FindNeighbors(pca_obj, dims = 1:pca_dims)

message(glue("Identifying clusters at resolution {cluster_res}..."))
# This groups them into numbered clusters (0, 1, 2, etc.)
cluster_obj <- FindClusters(cluster_obj, resolution = cluster_res)

# ========== 4. RUN UMAP ==========
message("Running UMAP for 2D visualization...")
umap_obj <- RunUMAP(cluster_obj, dims = 1:pca_dims)

# ========== 5. GENERATE PLOTS ==========
message("Generating UMAP plots...")

# Plot 1: Colored by Cluster (What cell types are there?)
p1 <- DimPlot(umap_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP by Cluster")

# Plot 2: Colored by Condition (Untreated vs Treated)
p2 <- DimPlot(umap_obj, reduction = "umap", group.by = "condition", pt.size = 0.5) +
  ggtitle(umap_title)

# Plot 2a: Plot them side-by-side instead of on top of each other!
p3 <- DimPlot(umap_obj, reduction = "umap", group.by = "condition", split.by = "condition")

# Combine them side-by-side using the patchwork library
combined_plot <- p1 + p2 + p3

pdf(path(plot_dir, "05_final_UMAP.pdf"), width = 12, height = 5)
print(combined_plot)
dev.off()

# ========== 6. SAVE FINAL OBJECT ==========
message("Saving final clustered UMAP Seurat object...")
saveRDS(umap_obj, path(obj_dir, "05_umap.rds"))

message("Step 05 Complete! Data clustered and UMAP generated.")

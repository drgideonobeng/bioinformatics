# scripts/06_run_dim_reduction.R
library(Seurat)
library(fs)
library(glue)

# 1. Pull dynamic variables from config
obj_dir     <- Sys.getenv("OBJ_DIR")
plot_dir    <- Sys.getenv("PLOT_DIR")
pca_dims    <- as.numeric(Sys.getenv("PCA_DIMS", unset = 15))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.8))

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the normalized object
load_path <- path(obj_dir, "05_seurat_normalized_scaled.rds")
message("Loading normalized & scaled data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Run PCA (Calculated only once!)
message("Running PCA...")
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)

# 3b. Generate the diagnostic Elbow Plot
message("Generating Elbow Plot to determine appropriate number of PCA dims")
dir_create(plot_dir)
pdf(path(plot_dir, "06_elbow_plot.pdf"), width = 8, height = 6)
print(ElbowPlot(seurat_obj, ndims = 50))
dev.off()
message("Elbow plot saved to: ", path(plot_dir, "06_elbow_plot.pdf"))

# 3c. The JackStraw Method (The Statistical Way to determine PCA dims)
# message("Generating JackStraw Plot to determine appropriate number of PCA dims")
# seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
# seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:20)
# pdf(path(plot_dir, "06_jackstraw_plot.pdf"), width = 8, height = 6) 
# print(JackStrawPlot(seurat_obj, dims = 1:20)) 
# dev.off()
# message("JackStraw plot saved to: ", path(plot_dir, "06_jackstraw_plot.pdf"))

# 4. Clustering (Finding the biological groups)
message(glue("Finding neighbors using top {pca_dims} PCA dimensions..."))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pca_dims)

message(glue("Finding clusters with resolution: {round(cluster_res, 2)}..."))
seurat_obj <- FindClusters(seurat_obj, resolution = cluster_res)

# 5. Run UMAP (For 2D Visualization)
message("Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pca_dims)

# 6. Save the fully processed object
save_path <- path(obj_dir, "06_seurat_clustered.rds")
saveRDS(seurat_obj, file = save_path)

message("Step 06 Complete. Clustered object saved to: ", save_path)

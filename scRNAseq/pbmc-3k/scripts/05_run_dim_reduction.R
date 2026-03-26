# scripts/05_run_dim_reduction.R
library(Seurat)
library(fs)

# 1. Pull dynamic variables from config
obj_dir <- Sys.getenv("OBJ_DIR")
pca_dims <- as.numeric(Sys.getenv("PCA_DIMS", unset = 10))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.5))

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the filtered and normalized object
load_path <- path(obj_dir, "02_pbmc_filtered_normalized.rds")
message("Loading filtered data from: ", load_path)
pbmc <- readRDS(load_path)

# 3. Run PCA (Principal Component Analysis)
message("Running PCA...")
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)

# 4. Clustering (Finding the biological groups)
# We use the top N dimensions defined in your config (usually 10)
message("Finding neighbors using ", pca_dims, " PCA dimensions...")
pbmc <- FindNeighbors(pbmc, dims = 1:pca_dims)

message("Finding clusters with resolution: ", cluster_res, "...")
pbmc <- FindClusters(pbmc, resolution = cluster_res)

# 5. Run UMAP (For 2D Visualization)
message("Running UMAP...")
pbmc <- RunUMAP(pbmc, dims = 1:pca_dims)

# 6. Save the fully processed object
save_path <- path(obj_dir, "03_pbmc_clustered.rds")
saveRDS(pbmc, file = save_path)

message("Success! Clustered object saved to: ", save_path)

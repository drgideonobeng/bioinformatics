#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(ggplot2)
library(fs)
library(glue)

# Fetch configurations
obj_dir   <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")
dim_start <- as.numeric(Sys.getenv("LSI_DIMS_START", unset = 2))
dim_end   <- as.numeric(Sys.getenv("LSI_DIMS_END", unset = 30))
clust_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.8))

message(glue("==============================================="))
message(glue(" Step 05: scATAC-seq LSI, UMAP & Clustering    "))
message(glue("==============================================="))

# 1. Load Filtered Object
read_path <- path(obj_dir, "04_atac_filtered.rds")
message(glue("Loading filtered object: {read_path}"))
seurat_obj <- readRDS(read_path)

# 2. Dimensionality Reduction (LSI)
message(glue("Performing TF-IDF Normalization..."))
seurat_obj <- RunTFIDF(seurat_obj)

message(glue("Finding Top Features..."))
seurat_obj <- FindTopFeatures(seurat_obj, min.cutoff = 'q0')

message(glue("Running Singular Value Decomposition (SVD)..."))
seurat_obj <- RunSVD(seurat_obj)

# 3. Non-linear Dimension Reduction and Clustering
lsi_dims <- dim_start:dim_end
message(glue("Running UMAP & Clustering using LSI Dimensions {dim_start} to {dim_end}..."))

seurat_obj <- RunUMAP(seurat_obj, reduction = 'lsi', dims = lsi_dims)
seurat_obj <- FindNeighbors(seurat_obj, reduction = 'lsi', dims = lsi_dims)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, algorithm = 3, resolution = clust_res)

# 4. Generate UMAP Plot
message(glue("Plotting UMAP..."))
p1 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.5) + 
      ggtitle("scATAC-seq PBMC Clustering") +
      theme_minimal()

plot_path <- path(plot_dir, "05_atac_umap.pdf")
ggsave(filename = plot_path, plot = p1, width = 8, height = 6)
message(glue("-> Saved UMAP plot to: {plot_path}"))

# 5. Save Object
save_path <- path(obj_dir, "05_atac_clustered.rds")
saveRDS(seurat_obj, file = save_path)
message("===============================================")
message(glue("-Step 05 Completed:Saved clustered object to: {save_path}"))
message("===============================================")

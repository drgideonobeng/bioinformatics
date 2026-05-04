#!/usr/bin/env Rscript

# scripts/02a_qc_and_integration.R

library(Seurat)
library(fs)
library(glue)
library(harmony)
library(tidyverse)
library(future)

# Increase threads/workers and access limit to 8GB
plan(multicore, workers = 8)
options(future.globals.maxSize = 8 * 1024^3)

# 1. Pull dynamic variables from config
matrix_dir  <- Sys.getenv("DATA_MATRIX_DIR")
obj_dir     <- Sys.getenv("OBJ_DIR")
plot_dir    <- Sys.getenv("PLOT_DIR")
mt_pattern  <- Sys.getenv("MT_PATTERN", unset = "^mt-")
min_cells   <- as.numeric(Sys.getenv("MIN_CELLS", unset = 3))
min_genes   <- as.numeric(Sys.getenv("MIN_GENES", unset = 200))
max_genes   <- as.numeric(Sys.getenv("MAX_GENES", unset = 6000))
max_mt      <- as.numeric(Sys.getenv("MAX_MT_PERCENT", unset = 5))
pca_dims    <- as.numeric(Sys.getenv("PCA_DIMS", unset = 20))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.5))

if (matrix_dir == "") stop("DATA_MATRIX_DIR is empty. Please run 'source config.sh'.")

# 2. Load the sample metadata
meta <- read_csv("sample_metadata.csv", show_col_types = FALSE)
seurat_list <- list() # initialize the list

message("--- Phase 1: Loading and Filtering Individual Samples ---")
# Phase 1: Loading and Filtering Individual Samples
message("--- Phase 1: Loading and Filtering Individual Samples ---")
for (i in 1:nrow(meta)) {
  # Extract values for this iteration
  sample_id      <- meta$sample_id[i]
  genotype <- meta$genotype[i]
  batch    <- meta$batch[i]
  
  sample_mat_dir <- path(matrix_dir, sample_id)
  message(glue("Processing {sample_id} ({genotype})..."))
  
  # Check if directory exists before loading
  if (!dir_exists(sample_mat_dir)) stop(glue("Directory not found: {sample_mat_dir}"))

  # Load and Create Object
  counts <- Read10X(data.dir = sample_mat_dir)
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = min_cells, min.features = min_genes)
  
  # Assign metadata
  obj$genotype <- genotype
  obj$batch    <- batch
  
  # QC and Filtering
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj <- subset(obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)
  
  # Add to our list
  seurat_list[[sample_id]] <- obj
}

# Phase 2: Merging Objects
message("\n--- Phase 2: Merging Objects ---")
# Merge the first object with the rest of the list
merged_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
message(glue("Giant merged object created with {ncol(merged_obj)} total cells."))

# Phase 3 : Standard Pre-processing (Normalize, Find variable features, Scale, Run PCA)
message("\n--- Phase 3: Standard Pre-processing ---")
merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000)
merged_obj <- ScaleData(merged_obj, vars.to.regress = "percent.mt")
merged_obj <- RunPCA(merged_obj, npcs = 50, verbose = FALSE)

message("\n--- Generating PCA Elbow Plot ---")
dir_create(plot_dir)
pdf(path(plot_dir, "02_elbow_plot.pdf"), width = 8, height = 6)
print(ElbowPlot(merged_obj, ndims = 50))
dev.off()
message("Elbow plot saved to plots directory.")
# ------------------------------------------------------------

# Phase 4: Harmony Integration
message("\n--- Phase 4: Harmony Integration ---")
# group.by.vars tells Harmony which metadata column contains the technical batches
merged_obj <- RunHarmony(merged_obj, group.by.vars = "batch", plot_convergence = FALSE)

# Phase 5: Dimensionality Reduction & Clustering
message("\n--- Phase 5: Dimensionality Reduction & Clustering ---")
# CRITICAL: We must specify reduction = "harmony" instead of "pca"
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:pca_dims)
merged_obj <- FindNeighbors(merged_obj, reduction = "harmony", dims = 1:pca_dims)
merged_obj <- FindClusters(merged_obj, resolution = cluster_res)

# Saving Visuals and Objects
message("\n--- Saving Visuals and Objects ---")
dir_create(plot_dir)
dir_create(obj_dir)

# Save UMAP split by Genotype to verify integration
pdf(path(plot_dir, "02_umap_split.pdf"), width = 12, height = 6)
p1 <- DimPlot(merged_obj, reduction = "umap", split.by = "genotype", group.by = "seurat_clusters", label = TRUE) + 
      ggtitle("Harmony Integration: Split by Genotype")
print(p1)
dev.off()

# Save final object
save_path <- path(obj_dir, "02_seurat_integrated.rds")
saveRDS(merged_obj, file = save_path)

message("\nStep 02a:QC and Harmony Integration complete. Visuals saved to plots directory, object saved to: ", save_path)

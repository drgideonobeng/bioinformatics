#!/usr/bin/env Rscript

# ========== 0. LOAD PACKAGES ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(fs)
  library(glue)
  library(harmony)
  library(tidyverse)
  library(future)
})

# ========== 1. CONFIGURATION & ENVIRONMENT ==========
matrix_dir  <- Sys.getenv("DATA_MATRIX_DIR")
obj_dir     <- Sys.getenv("OBJ_DIR")
plot_dir    <- Sys.getenv("PLOT_DIR")
mt_pattern  <- Sys.getenv("MT_PATTERN", unset = "^MT-")
min_cells   <- as.numeric(Sys.getenv("MIN_CELLS", unset = 3))
min_genes   <- as.numeric(Sys.getenv("MIN_GENES", unset = 200))
max_genes   <- as.numeric(Sys.getenv("MAX_GENES", unset = 10000))
max_mt      <- as.numeric(Sys.getenv("MAX_MT_PERCENT", unset = 10))
pca_dims    <- as.numeric(Sys.getenv("PCA_DIMS", unset = 20))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.5))

if (matrix_dir == "") stop("DATA_MATRIX_DIR is empty. Please run 'source config.sh'.")

# Setup parallel processing
n_cores <- as.numeric(Sys.getenv("THREADS", unset = 8))
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 8000 * 1024^2)

# ========== 2. LOAD SAMPLE METADATA ==========
meta <- read_csv("sample_metadata.csv", show_col_types = FALSE)
seurat_list <- list() 

# ========== 3. PHASE 1: LOAD & FILTER SAMPLES ==========
message("--- Phase 1: Loading and Filtering Individual Samples ---")
for (i in 1:nrow(meta)) {
  sid      <- meta$sample_id[i]
  sample_type <- meta$sample_type[i]
  batch    <- meta$batch[i]
  
  sample_mat_dir <- path(matrix_dir, sid)
  message(glue("Processing {sid} ({sample_type})..."))
  
  if (!dir_exists(sample_mat_dir)) stop(glue("Directory not found: {sample_mat_dir}"))

  counts <- Read10X(data.dir = sample_mat_dir)
  obj <- CreateSeuratObject(counts = counts, project = sid, min.cells = min_cells, min.features = min_genes)
  
  obj$sample_type <- sample_type
  obj$batch    <- batch
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj <- subset(obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)
  
  seurat_list[[sid]] <- obj
}

# ========== 4. PHASE 2: MERGE & FACTORIZE ==========
message("\n--- Phase 2: Merging Objects ---")
merged_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
message(glue("Merged object created with {ncol(merged_obj)} total cells."))

# Convert metadata variables to factors
merged_obj$batch <- as.factor(merged_obj$batch)
merged_obj$sample_type <- as.factor(merged_obj$sample_type)

# Count unique batches to determine if integration is mathematically possible
n_batches <- length(unique(merged_obj$batch))
message(glue("Detected {n_batches} unique batch(es) in metadata."))

# ========== 5. PHASE 3: PRE-PROCESSING ==========
message("\n--- Phase 3: Standard Pre-processing ---")
merged_obj <- NormalizeData(merged_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

message(glue("Scaling data and regressing out percent.mt using {n_cores} cores..."))
merged_obj <- ScaleData(merged_obj, vars.to.regress = "percent.mt", verbose = TRUE)
merged_obj <- RunPCA(merged_obj, npcs = 50, verbose = FALSE)

message("\n--- Generating PCA Elbow Plot ---")
dir_create(plot_dir)
pdf(path(plot_dir, "02_elbow_plot.pdf"), width = 8, height = 6)
print(ElbowPlot(merged_obj, ndims = 50))
invisible(dev.off())

# ========== 6. PHASE 4: CONDITIONAL HARMONY INTEGRATION ==========
if (n_batches > 1) {
  message("\n--- Phase 4: Harmony Integration ---")
  merged_obj <- RunHarmony(merged_obj, group.by.vars = "batch", plot_convergence = FALSE)
  reduction_to_use <- "harmony"
} else {
  message("\n--- Skipping Phase 4: Only one batch detected. Defaulting to PCA for UMAP. ---")
  reduction_to_use <- "pca"
}

# ========== 7. PHASE 5: DIMENSIONALITY REDUCTION & CLUSTERING ==========
message("\n--- Phase 5: Dimensionality Reduction & Clustering ---")
# Uses either "harmony" or "pca" dynamically based on the Phase 4 check
merged_obj <- RunUMAP(merged_obj, reduction = reduction_to_use, dims = 1:pca_dims, verbose = FALSE)
merged_obj <- FindNeighbors(merged_obj, reduction = reduction_to_use, dims = 1:pca_dims, verbose = FALSE)
merged_obj <- FindClusters(merged_obj, resolution = cluster_res, verbose = FALSE)

# ========== 8. SAVE VISUALS & OBJECT ==========
message("\n--- Saving Visuals and Objects ---")
dir_create(obj_dir)

# Save UMAP split by Genotype
pdf(path(plot_dir, "02_integrated_umap_split.pdf"), width = 12, height = 6)
p1 <- DimPlot(merged_obj, reduction = "umap", split.by = "sample_type", group.by = "seurat_clusters", label = TRUE) + 
      ggtitle(glue("UMAP Integration: Split by Genotype (Reduction: {reduction_to_use})"))
print(p1)
invisible(dev.off())

save_path <- path(obj_dir, "02_seurat_harmony_integrated.rds")
saveRDS(merged_obj, file = save_path)

message("\nStep 02 Complete. Visuals saved to plots directory. Object saved to: ", save_path)

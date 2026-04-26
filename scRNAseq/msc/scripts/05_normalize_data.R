#!/usr/bin/env Rscript

# suppressPackageStartupMessages keeps your terminal clean from loading spam
suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(fs)
  library(glue)
})

# ========== 1. CONFIGURATION & ENVIRONMENT ==========
obj_dir <- Sys.getenv("OBJ_DIR")
if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# Setup parallel processing using Mac-friendly "multicore"
n_cores <- as.numeric(Sys.getenv("THREADS", unset = 8))
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 8000 * 1024^2) # 8GB limit per thread

message(glue("Starting Step 05: Normalization & Scaling ({n_cores} threads allocated)"))

# ========== 2. LOAD DATA ==========
load_path <- path(obj_dir, "04_seurat_filtered.rds")
message(glue("Loading filtered object from: {load_path}"))
seurat_obj <- readRDS(load_path)

# ========== 3. NORMALIZE & FIND VARIABLE FEATURES ==========
message("Normalizing data (LogNormalize)...")
seurat_obj <- NormalizeData(
  seurat_obj, 
  normalization.method = "LogNormalize", 
  scale.factor = 10000,
  verbose = FALSE # Silenced to keep terminal clean
)

message("Finding top 2,000 highly variable features...")
seurat_obj <- FindVariableFeatures(
  seurat_obj, 
  selection.method = "vst", 
  nfeatures = 2000,
  verbose = FALSE
)

# Extract and display the top 20 genes
top20_genes <- head(VariableFeatures(seurat_obj), 20)
message("Top 20 variable genes driving your dataset:")
print(top20_genes)

# ========== 4. SCALE DATA & REGRESS MT ==========
# This is the heavy step where parallelization kicks in
message(glue("Scaling all genes and regressing out 'percent.mt' using {n_cores} cores..."))

all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(
  seurat_obj, 
  features = all_genes, 
  vars.to.regress = "percent.mt",
  verbose = TRUE # Kept TRUE so you can see the progress bar
)

# ========== 5. SAVE PROGRESS ==========
save_path <- path(obj_dir, "05_seurat_normalized_scaled.rds")
message(glue("Saving normalized and scaled object to: {save_path}"))
saveRDS(seurat_obj, file = save_path)

message("Step 05 Complete!")

# scripts/03_normalize.R
library(Seurat)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 03...")
obj_dir <- Sys.getenv("OBJ_DIR")

# ========== 2. LOAD FILTERED OBJECT ==========
message("Loading QC-filtered dataset...")
filtered_obj <- readRDS(path(obj_dir, "02_qc_filtered.rds"))

# ========== 3. NORMALIZE DATA ==========
message("Normalizing RNA counts...")
# LogNormalize: Divides by total counts, multiplies by 10,000, and natural-log transforms
norm_obj <- NormalizeData(filtered_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# ========== 4. FIND VARIABLE FEATURES ==========
message("Identifying highly variable genes...")
# We look for the top 2000 genes that change the most between cells
norm_obj <- FindVariableFeatures(norm_obj, selection.method = "vst", nfeatures = 2000)

# ========== 5. SAVE NORMALIZED OBJECT ==========
message("Saving normalized Seurat object...")
saveRDS(norm_obj, path(obj_dir, "03_normalized.rds"))

message("Step 03 Complete! Data normalized and highly variable features identified.")

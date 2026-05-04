# scripts/05_normalize_data.R
library(Seurat)
library(fs)
library(future)

# Increase threads/workers to 8
plan(multicore, workers = 8)
# Increase global access limit per worker to 8GB
options(future.globals.maxSize= 8 * 1024^3)

message("\n==========Parallel processing activated==========")

# Pull dynamic variables from config.sh
obj_dir <- Sys.getenv("OBJ_DIR")
if (obj_dir == "") stop("OBJ_DIR is empty. Run source config.sh or via run_pipeline.sh")

# 1. Load the filtered object
message("\n=====loading seurat object=====")
load_path <- path(obj_dir, "04_seurat_filtered.rds")
message("Loading filtered object from: ", load_path)
seurat_obj <- readRDS(load_path)

# 2. Normalize the data
message("Normalizing data (LogNormalize)...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 3. Identify highly variable features
message("Finding highly variable features...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 4. Scale the data
message("Scaling data...")
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes, vars.to.regress = "percent.mt")

# 5. Save the normalized and scaled object
save_path <- path(obj_dir, "05_seurat_normalized_scaled.rds")
saveRDS(seurat_obj, file = save_path)

message("Step 05 Complete. Normalized & scaled object saved to: ", save_path)

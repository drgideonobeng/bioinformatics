# scripts/05_normalize_data.R
library(Seurat)
library(fs)

obj_dir <- Sys.getenv("OBJ_DIR")
if (obj_dir == "") stop("OBJ_DIR is empty. Run source config.sh or via run_pipeline.sh")

# 1. Load the filtered object
load_path <- path(obj_dir, "02_seurat_filtered.rds")
message("Loading filtered object from: ", load_path)
seurat_obj <- readRDS(load_path)

# 2. Normalize the data
message("Normalizing data (LogNormalize)...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

# 3. Identify highly variable features
message("Finding highly variable features...")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# 3b. Identify the 10 most highly varible genes
message("Top 10 variable genes")
top10 <- head(VariableFeatures(seurat_obj), 10)
print(top10)

# 4. Scale the data
message("Scaling data...")
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes, vars.to.regress = "percent.mt")

# 5. Save the normalized and scaled object
save_path <- path(obj_dir, "03_seurat_normalized.rds")
saveRDS(seurat_obj, file = save_path)

message("Step 05 Complete. Normalized object saved to: ", save_path)

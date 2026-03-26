# scripts/04_filter_cells.R
library(Seurat)
library(fs)

# 1. Pull dynamic variables from config (with fallbacks)
obj_dir   <- Sys.getenv("OBJ_DIR")
min_genes <- as.numeric(Sys.getenv("MIN_GENES", unset = 200))
max_genes <- as.numeric(Sys.getenv("MAX_GENES", unset = 2500))
max_mt    <- as.numeric(Sys.getenv("MAX_MT_PERCENT", unset = 5))

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the unfiltered object
load_path <- path(obj_dir, "01_pbmc_unfiltered.rds")
message("Loading unfiltered object...")
pbmc <- readRDS(load_path)

# Record starting cell count for our logs
pre_filter_cells <- ncol(pbmc)

# 3. Filter the cells
message(sprintf("Filtering: Genes > %d & Genes < %d & Mito < %d%%", min_genes, max_genes, max_mt))

# The actual subset command
pbmc <- subset(pbmc, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)

# Record ending cell count
post_filter_cells <- ncol(pbmc)
message("Retained ", post_filter_cells, " out of ", pre_filter_cells, " cells.")

# 4. Normalize the data
message("Normalizing data (LogNormalize)...")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# 5. Identify highly variable features (Genes that drive biological differences)
message("Finding highly variable features...")
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# 6. Scale the data (Mean = 0, Variance = 1 for PCA)
message("Scaling data...")
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 7. Save the processed object
save_path <- path(obj_dir, "02_pbmc_filtered_normalized.rds")
saveRDS(pbmc, file = save_path)

message("Success! Filtered and normalized object saved to: ", save_path)

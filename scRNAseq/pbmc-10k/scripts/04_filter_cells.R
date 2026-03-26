# scripts/04_filter_cells.R
library(Seurat)
library(fs)
library(glue)

# 1. Pull dynamic variables from config
obj_dir   <- Sys.getenv("OBJ_DIR")
min_genes <- as.numeric(Sys.getenv("MIN_GENES", unset = 200))
max_genes <- as.numeric(Sys.getenv("MAX_GENES", unset = 4000))
max_mt    <- as.numeric(Sys.getenv("MAX_MT_PERCENT", unset = 5))

if (obj_dir == "") stop("OBJ_DIR is empty. Run source config.sh first or  via run_pipeline.sh")

# 2. Load the unfiltered object
load_path <- path(obj_dir, "01_seurat_unfiltered.rds")
message("Loading unfiltered object from: ", load_path)
seurat_obj <- readRDS(load_path)

pre_filter_cells <- ncol(seurat_obj)

# 3. Filter the cells
message(glue("Filtering: Genes > {min_genes} & Genes < {max_genes} & Mito < {max_mt}%"))
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)

post_filter_cells <- ncol(seurat_obj)
message("Retained ", post_filter_cells, " out of ", pre_filter_cells, " cells.")

# 4. Save the cleanly filtered (but un-normalized) object
save_path <- path(obj_dir, "02_seurat_filtered.rds")
saveRDS(seurat_obj, file = save_path)

message("Step 04 Complete. Filtered object saved to: ", save_path)

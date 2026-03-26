# scripts/01_load_and_merge_data.R
library(Seurat)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 01...")
data_dir_untreated <- Sys.getenv("DATA_DIR_UNTREATED")
data_dir_treated   <- Sys.getenv("DATA_DIR_TREATED")
obj_dir          <- Sys.getenv("OBJ_DIR")
project_name     <- Sys.getenv("PROJECT_NAME")

min_cells        <- as.numeric(Sys.getenv("MIN_CELLS"))
min_features     <- as.numeric(Sys.getenv("MIN_FEATURES_BASE"))

# ========== 2. LOAD & MERGE DATA ==========
message("Reading raw 10x matrices...")
untreated_matrix <- Read10X(data.dir = data_dir_untreated)
treated_matrix   <- Read10X(data.dir = data_dir_treated)

message("Creating base Seurat objects...")
untreated_seurat <- CreateSeuratObject(counts = untreated_matrix, project = "Untreated", min.cells = min_cells, min.features = min_features)
treated_seurat   <- CreateSeuratObject(counts = treated_matrix, project = "Treated", min.cells = min_cells, min.features = min_features)

message(glue("{ncol(untreated_seurat)} & {ncol(treated_seurat)} cells counted in the untreated and treated group respectively "))

# Add metadata 
untreated_seurat$condition <- "Untreated"
treated_seurat$condition <- "Treated"

message("Merging datasets...")
merged_obj <- merge(untreated_seurat, y = treated_seurat, add.cell.ids = c("H", "T"), project = project_name)

# ========== 3. SAVE RAW OBJECT ==========
message("Saving raw merged Seurat object...")
saveRDS(merged_obj, path(obj_dir, "01_raw_merged.rds"))

message(glue("Step 01 Complete! Raw dataset merged with {ncol(merged_obj)} total cells."))

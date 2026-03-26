# scripts/02_create_seurat_obj.R

library(Seurat)
library(fs)

# 1. Pull dynamic variables from config (with a fallback default)
data_dir <- Sys.getenv("DATA_DIR")
obj_dir  <- Sys.getenv("OBJ_DIR")
min_genes <- as.numeric(Sys.getenv("MIN_GENES", unset = 200)) 

# Safety check
if (data_dir == "") stop("DATA_DIR is empty. Please run 'source config.sh' first.")

message("Loading 10x Data from: ", data_dir)
message("Applying MIN_GENES threshold: ", min_genes)

# 2. Load data and create the object
pbmc.data <- Read10X(data.dir = data_dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = min_genes)

# 3. Calculate Mitochondrial Percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 4. Save the object cleanly
dir_create(obj_dir) # fs handles this silently if the folder already exists
save_path <- path(obj_dir, "01_pbmc_unfiltered.rds")
saveRDS(pbmc, file = save_path)

message("Success! Unfiltered Seurat object saved to: ", save_path)

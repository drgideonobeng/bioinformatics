# scripts/02_create_seurat_obj.R

library(Seurat)
library(fs)

# 1. Pull dynamic variables from config
matrix_dir <- Sys.getenv("DATA_MATRIX_DIR")
project    <- Sys.getenv("PROJECT_NAME")
mt_pattern <- Sys.getenv("MT_PATTERN")
obj_dir    <- Sys.getenv("OBJ_DIR")

# Dynamically pull the prefilter numbers
min_cells    <- as.numeric(Sys.getenv("MIN_CELLS"))
min_features <- as.numeric(Sys.getenv("MIN_FEATURES_BASE"))

# Safety check
if (matrix_dir == "") stop("DATA_MATRIX_DIR is empty. Please run 'source config.sh'.")

# 2. Load data
message("Loading 10x matrix from: ", matrix_dir)
raw.data <- Read10X(data.dir = matrix_dir)

# 2b. view n[rows, columns]Output will be [Number of Genes][Number of Cells]
n_genes <- nrow(raw.data)
n_cells <- ncol(raw.data)
message("Matrix contains ",n_genes, " Genes & ",  n_cells, " Cells" )

# 3. Create Seurat object 
message("Creating Seurat Object for", project, "...")
seurat_obj <- CreateSeuratObject(counts = raw.data, 
                           project = project, 
                           min.cells = min_cells, 
                           min.features = min_features)

# 4. Calculate Mitochondrial Percentage
message("Calculating MT% using pattern: ", mt_pattern)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)

# 4b. Show QC metrics stored in object metadata
message("Showing QC metrics for the first 20 cells")
head(seurat_obj@meta.data, 20)

# 5. Save the object cleanly
dir_create(obj_dir) # fs handles this silently if the folder already exists
save_path <- path(obj_dir, "01_seurat_unfiltered.rds")
saveRDS(seurat_obj, file = save_path)

message("Success! Unfiltered Seurat object saved to: ", save_path)

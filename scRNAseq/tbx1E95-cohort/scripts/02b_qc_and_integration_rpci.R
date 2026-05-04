#!/usr/bin/env Rscript

# scripts/02b_qc_and_integration_rpci.R

library(Seurat)
library(fs)
library(glue)
library(tidyverse)
library(future)

# Load RISC (Install from GitHub if missing: remotes::install_github("bioinfoDZ/RISC"))
if (!requireNamespace("RISC", quietly = TRUE)) {
  stop("The RISC package is required for RPCI integration. Please install it via GitHub.")
}
library(RISC)

# Increase threads/workers and access limit to 8GB
plan(multicore, workers = 8)
options(future.globals.maxSize = 8 * 1024^3)
message("\n==========Parallel processing activated==========\n")

# 1. Pull dynamic variables from config
matrix_dir  <- Sys.getenv("DATA_MATRIX_DIR")
obj_dir     <- Sys.getenv("OBJ_DIR")
plot_dir    <- Sys.getenv("PLOT_DIR")
mt_pattern  <- Sys.getenv("MT_PATTERN", unset = "^mt-")
min_cells   <- as.numeric(Sys.getenv("MIN_CELLS", unset = 3))
min_genes   <- as.numeric(Sys.getenv("MIN_GENES", unset = 200))
max_genes   <- as.numeric(Sys.getenv("MAX_GENES", unset = 6000))
max_mt      <- as.numeric(Sys.getenv("MAX_MT_PERCENT", unset = 5))
pca_dims    <- as.numeric(Sys.getenv("PCA_DIMS", unset = 20))
cluster_res <- as.numeric(Sys.getenv("CLUSTER_RES", unset = 0.5))

if (matrix_dir == "") stop("DATA_MATRIX_DIR is empty. Please run 'source config.sh'.")

# 2. Load the sample metadata
meta <- read_csv("sample_metadata.csv", show_col_types = FALSE)
seurat_list <- list()
risc_list   <- list()

message("--- Phase 1: Loading, QC, and RISC Preparation ---")
for (i in 1:nrow(meta)) {
  sample_id  <- meta$sample_id[i]
  genotype   <- meta$genotype[i]
  batch      <- meta$batch[i]
  
  sample_mat_dir <- path(matrix_dir, sample_id)
  message(glue("Processing {sample_id} ({genotype})..."))
  
  if (!dir_exists(sample_mat_dir)) stop(glue("Directory not found: {sample_mat_dir}"))

  # A. Load and Create Seurat Object (for consistent QC)
  counts <- Read10X(data.dir = sample_mat_dir)
  obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = min_cells, min.features = min_genes)
  
  obj$genotype <- genotype
  obj$batch    <- batch
  
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj <- subset(obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)
  
  # Add cell prefixes so names don't collide when merging
  obj <- RenameCells(obj, add.cell.id = sample_id)
  
  # B. Convert to RISC format
  message(glue("Converting {sample_id} to RISC format..."))
  raw_counts <- GetAssayData(obj, layer = "counts")
  cell_meta  <- obj@meta.data
  
  # Initialize RISC object, Normalize, and find Dispersions (Variable Genes)
  r_obj <- readsc(raw_counts, cell = cell_meta, gene = data.frame(Symbol = rownames(raw_counts), row.names = rownames(raw_counts)))
  r_obj <- scFilter(r_obj)
  r_obj <- scNormalize(r_obj)
  r_obj <- scDisperse(r_obj)
  
  # CRITICAL SYNC: RISC dropped some cells during scFilter. 
  # Force Seurat to drop those exact same cells so they are 1:1.
  obj <- subset(obj, cells = rownames(r_obj@coldata))
  
  seurat_list[[sample_id]] <- obj
  risc_list[[sample_id]]   <- r_obj
}

message("\n--- Phase 2: RPCI Integration via RISC ---")
# Identify which dataset is our Reference.
ref_index <- which(meta$genotype == "WT" | meta$genotype == "Wild-Type")[1]
if (is.na(ref_index)) ref_index <- 1 # Fallback just in case

message(glue("Using Sample '{meta$sample_id[ref_index]}' as the Reference Dataset for RPCI..."))

# Reorder BOTH lists so the reference (WT) is the very first object.
# This ensures Seurat merges in the exact same order that RISC integrates.
risc_list   <- c(risc_list[ref_index], risc_list[-ref_index])
seurat_list <- c(seurat_list[ref_index], seurat_list[-ref_index])

# Run the asymmetrical integration
integrated_risc <- scMultiIntegrate(
  objects = risc_list, 
  align = "Predict", 
  npc = 50
)

message("\n--- Phase 3: Porting RPCI Embeddings back to Seurat ---")
# Merge the original Seurat objects in the newly reordered order
merged_seurat <- merge(x = seurat_list[[1]], y = seurat_list[-1])
merged_seurat <- NormalizeData(merged_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat <- FindVariableFeatures(merged_seurat)

# Extract the corrected RPCI coordinates from the RISC object
rpci_embeddings <- integrated_risc@DimReduction$cell.pls

# Because we synced the filtering in Phase 1, and the ordering in Phase 2,
# the cells align perfectly. We can safely stamp the names back onto the matrix.
rownames(rpci_embeddings) <- colnames(merged_seurat)

# Create a custom Dimensional Reduction slot in Seurat named "rpci"
merged_seurat[["rpci"]] <- CreateDimReducObject(
  embeddings = as.matrix(rpci_embeddings), 
  key = "RPCI_", 
  assay = DefaultAssay(merged_seurat)
)

message("\n--- Phase 4: Dimensionality Reduction & Clustering (RPCI) ---")
# CRITICAL: We instruct Seurat to build the UMAP/Clusters using the new "rpci" matrix
merged_seurat <- RunUMAP(merged_seurat, reduction = "rpci", dims = 1:pca_dims)
merged_seurat <- FindNeighbors(merged_seurat, reduction = "rpci", dims = 1:pca_dims)
merged_seurat <- FindClusters(merged_seurat, resolution = cluster_res)

message("\n--- Saving Visuals and Objects ---")
dir_create(plot_dir)
dir_create(obj_dir)

# Save UMAP split by Genotype to verify integration
pdf(path(plot_dir, "02b_rpci_umap_split.pdf"), width = 12, height = 6)
p1 <- Seurat::DimPlot(merged_seurat, reduction = "umap", split.by = "genotype", group.by = "seurat_clusters", label = TRUE) + 
      ggtitle("RPCI Integration: Reference Preserved Map")
print(p1)
invisible(dev.off())

# Save final object
save_path <- path(obj_dir, "02_seurat_integrated.rds")
saveRDS(merged_seurat, file = save_path)

message("\nStep 02b: QC and RPCI Integration complete. Visuals saved to plots directory, object saved to: ", save_path)

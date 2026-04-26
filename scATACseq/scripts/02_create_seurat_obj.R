#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(fs)
library(glue)

# Fetch configurations with fallbacks
matrix_dir   <- Sys.getenv("DATA_MATRIX_DIR", unset = "data/raw_matrix")
project      <- Sys.getenv("PROJECT_NAME", unset = "scATAC_Project")
obj_dir      <- Sys.getenv("OBJ_DIR", unset = "results/objects")
min_cells    <- as.numeric(Sys.getenv("MIN_CELLS", unset = 10))
min_features <- as.numeric(Sys.getenv("MIN_FEATURES", unset = 200))
genome_build <- Sys.getenv("GENOME_BUILD", unset = "hg38")

h5_path   <- path(matrix_dir, "matrix.h5")
frag_path <- path(matrix_dir, "fragments.tsv.gz")

message(glue("Loading ATAC matrix from {h5_path}..."))
counts <- Read10X_h5(filename = h5_path)

message(glue("Building ChromatinAssay..."))
chrom_assay <- CreateChromatinAssay(
  counts       = counts,
  sep          = c(":", "-"),
  fragments    = frag_path,
  min.cells    = min_cells,
  min.features = min_features
)

seurat_obj <- CreateSeuratObject(
  counts  = chrom_assay,
  assay   = "peaks",
  project = project
)

message(glue("Extracting gene annotations and setting genome to {genome_build}..."))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- genome_build

Annotation(seurat_obj) <- annotations

dir_create(obj_dir) 
save_path <- path(obj_dir, "02_atac_unfiltered.rds")
saveRDS(seurat_obj, file = save_path)
message(glue("Success! Unfiltered ATAC object saved to: {save_path}"))

#!/usr/bin/env Rscript

# ========== 0. LOAD PACKAGES ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(fs)
  library(dplyr)
  library(readr) # Modern, faster CSV writing
  library(glue)
  library(future)
})

# ========== 1. CONFIGURATION & ENVIRONMENT ==========
obj_dir   <- Sys.getenv("OBJ_DIR")
table_dir <- Sys.getenv("TABLE_DIR")

if (obj_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# Setup parallel processing for FindAllMarkers
n_cores <- as.numeric(Sys.getenv("THREADS", unset = 8))
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 8000 * 1024^2)

# ========== 2. LOAD DATA ==========
load_path <- path(obj_dir, "02_seurat_harmony_integrated.rds")
message(glue("Loading integrated object from: {load_path}"))
seurat_obj <- readRDS(load_path)

# ========== 3. SEURAT V5 COMPATIBILITY ==========
message("Joining RNA layers for Seurat V5 compatibility...")
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- JoinLayers(seurat_obj)

# ========== 4. CALCULATE MARKERS ==========
message(glue("Calculating marker genes for all clusters across {n_cores} cores..."))
message("(This might take a few minutes depending on cell count)")

markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  verbose = TRUE
)

# ========== 5. PROCESS & EXPORT RESULTS ==========
dir_create(table_dir)
save_path <- path(table_dir, "03_harmony_cluster_markers.csv")
write_csv(markers, file = save_path)

# Extract top 5 genes per cluster safely
message("\n--- TOP 10 MARKER GENES PER CLUSTER ---")
top_genes <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>%
  ungroup() # Clean up grouping structure

top_genes_path <- path(table_dir, "03_top_markers_cheatsheet.csv")
write_csv(top_genes, file = top_genes_path)

# Print a clean sneak peek using dplyr
print(top_genes %>% select(cluster, gene, avg_log2FC), n = Inf)

message(glue("\nStep 03 Complete! Full marker list saved to: {save_path}"))

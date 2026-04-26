#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(fs)
  library(dplyr)
  library(readr) # Modern, faster data writing
  library(glue)
})

# ========== 1. CONFIGURATION & ENVIRONMENT ==========
obj_dir <- Sys.getenv("OBJ_DIR")
if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# Setup parallel processing for FindAllMarkers
n_cores <- as.numeric(Sys.getenv("THREADS", unset = 8))
plan("multicore", workers = n_cores)
options(future.globals.maxSize = 8000 * 1024^2)

message(glue("Starting Step 08: Marker Identification ({n_cores} threads allocated)"))

# ========== 2. LOAD DATA ==========
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message(glue("Loading clustered object from: {load_path}"))
seurat_obj <- readRDS(load_path)

# ========== 3. CALCULATE MARKERS ==========
message("Calculating marker genes for all clusters...")
message("(Distributing differential expression testing across your CPU cores)")

# Find markers across all clusters simultaneously
markers <- FindAllMarkers(
  seurat_obj, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  verbose = TRUE
)

# ========== 4. PROCESS & EXPORT RESULTS ==========
# Ensure the tables directory exists
tables_dir <- path(Sys.getenv("ROOT_DIR", unset = "."), "results", "tables")
dir_create(tables_dir)

# Save the full dataset
save_path <- path(tables_dir, "08_cluster_markers.csv")
write_csv(markers, file = save_path)

# Generate a top 10 cheat sheet
message("\n--- TOP 10 MARKER GENES PER CLUSTER ---")
top_genes <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  ungroup() # Clean up grouping structure

top_genes_path <- path(tables_dir, "08_top_cluster_markers_cheatsheet.csv")
write_csv(top_genes, file = top_genes_path)
 
# Print sneak peek using dplyr select
print(top_genes %>% select(cluster, gene, avg_log2FC), n = Inf)

message(glue("\nStep 08 Complete!"))
message(glue("Full marker list saved to: {save_path}"))
message(glue("Top genes cheat sheet saved to: {top_genes_path}"))

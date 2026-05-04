#!/usr/bin/env Rscript

# scripts/03_find_all_markers.R

library(Seurat)
library(fs)
library(dplyr)
library(glue)
library(future)

# Increase thread to 8 & access limit to 8GB
plan(multicore, workers = 8)
options(future.globals.maxSize = 8 * 1024^3)

# 1. Pull directories from config
obj_dir   <- Sys.getenv("OBJ_DIR")
table_dir <- Sys.getenv("TABLE_DIR")

if (obj_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the Harmony-integrated object
load_path <- path(obj_dir, "02_seurat_integrated.rds")
message("Loading integrated object from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. CRITICAL FOR SEURAT V5: Join the RNA layers before differential expression
message("Joining RNA layers for Seurat V5 compatibility...")
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- JoinLayers(seurat_obj)

# 4. Find all markers
message("Calculating marker genes for all clusters...")

# only.pos = TRUE means we only care about genes that are turned UP in a cluster
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

# 5. Save the results to a CSV file
dir_create(table_dir)
save_path <- path(table_dir, "03_cluster_markers.csv")
write.csv(markers, file = save_path, row.names = FALSE)

# 6. Print a sneak peek to the terminal and save a cheatsheet
message("\n--- TOP 10 MARKER GENES PER CLUSTER ---")
top_genes <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 10, order_by = avg_log2FC)

top_genes_path <- path(table_dir, "03_top_markers.csv")
write.csv(top_genes, file = top_genes_path, row.names = FALSE)
 
print(top_genes[, c("cluster", "gene", "avg_log2FC")], n= Inf)

message("\nStep 03 Complete! Full marker list saved to: ", save_path)

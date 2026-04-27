#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(fs)
library(glue)
library(readr) # tidyverse package for fast, clean CSV writing
library(future)

# Parrallel processing
plan("multicore", workers = 8)
options(future.globals.maxSize = 8000 * 1024^2)

message("Parallel processing enabled via future (M3 Max Optimized!)")

# Fetch configurations
obj_dir    <- Sys.getenv("OBJ_DIR", unset = "results/objects")
tables_dir <- Sys.getenv("TABLES_DIR", unset = "results/tables")

dir_create(tables_dir)

message(glue("==============================================="))
message(glue(" Step 07: Finding Marker Peaks (DA Analysis)   "))
message(glue("==============================================="))

# 1. Load the fully annotated object
read_path <- path(obj_dir, "06_atac_annotated.rds")
message(glue("Loading annotated ATAC object: {read_path}"))
atac_obj <- readRDS(read_path)

# 2. Set the active identities to our newly transferred RNA labels
Idents(atac_obj) <- "predicted.id"

# 3. Find Differentially Accessible (DA) Peaks
message(glue("Calculating Differentially Accessible peaks across all cell types..."))
message(glue("(This can take several minutes due to the massive size of the genome)"))

# Gold Standard ATAC Practice: LR test regressing out depth
da_peaks <- FindAllMarkers(
  object = atac_obj,
  only.pos = TRUE,         # We only care about regions that are OPENING, not closing
  min.pct = 0.05,          # The peak must be present in at least 5% of the target cells
  test.use = 'LR',         # Logistic Regression
  latent.vars = 'nCount_peaks' # Regress out sequencing depth bias
)

# 4. Save the results
save_path <- path(tables_dir, "07_marker_peaks.csv")
write_csv(da_peaks, save_path)
message(glue("-> Saved Marker Peaks table to: {save_path}"))
message(glue("==============================================="))

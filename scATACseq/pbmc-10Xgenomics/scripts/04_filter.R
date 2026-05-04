#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(fs)
library(glue)

obj_dir <- Sys.getenv("OBJ_DIR", unset = "results/objects")

# Dynamically fetch filtering thresholds from config.sh
min_peaks     <- as.numeric(Sys.getenv("MIN_PEAKS", unset = 3000))
max_peaks     <- as.numeric(Sys.getenv("MAX_PEAKS", unset = 100000))
min_tss       <- as.numeric(Sys.getenv("MIN_TSS_ENRICHMENT", unset = 2))
max_nucleo    <- as.numeric(Sys.getenv("MAX_NUCLEOSOME_SIGNAL", unset = 4))

message(glue("==============================================="))
message(glue(" Step 04: Filtering scATAC-seq Data            "))
message(glue("==============================================="))

# 1. Load object
read_path <- path(obj_dir, "03_atac_qc_calculated.rds")
message(glue("Loading object: {read_path}"))
seurat_obj <- readRDS(read_path)

pre_filter_cells <- ncol(seurat_obj)
message(glue("Pre-filter cell count: {pre_filter_cells}"))

# 2. Apply dynamic subset filters
message(glue("Applying subset filters:"))
message(glue(" - Peaks: {min_peaks} to {max_peaks}"))
message(glue(" - Min TSS: {min_tss}"))
message(glue(" - Max Nucleosome Signal: {max_nucleo}"))

seurat_filtered <- subset(
  x = seurat_obj,
  subset = nCount_peaks > min_peaks &
           nCount_peaks < max_peaks &
           TSS.enrichment > min_tss &
           nucleosome_signal < max_nucleo
)

post_filter_cells <- ncol(seurat_filtered)
message(glue("Post-filter cell count: {post_filter_cells}"))
message(glue("Total cells removed: {pre_filter_cells - post_filter_cells}"))

# 3. Save clean object
save_path <- path(obj_dir, "04_atac_filtered.rds")
saveRDS(seurat_filtered, file = save_path)
message(glue("-> Saved filtered object to: {save_path}"))
message(glue("==============================================="))

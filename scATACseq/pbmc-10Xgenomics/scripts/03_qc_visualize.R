#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(ggplot2)
library(fs)
library(glue)

obj_dir  <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")

dir_create(plot_dir)

message(glue("==============================================="))
message(glue(" Step 03: scATAC-seq Quality Control Metrics   "))
message(glue("==============================================="))

# 1. Load object
read_path <- path(obj_dir, "02_atac_unfiltered.rds")
message(glue("Loading unfiltered ATAC object from {read_path}..."))
seurat_obj <- readRDS(read_path)

# 2. Calculate Metrics
message(glue("Calculating Nucleosome Signal..."))
seurat_obj <- NucleosomeSignal(object = seurat_obj)

message(glue("Calculating TSS Enrichment Score (This may take a few minutes)..."))
seurat_obj <- TSSEnrichment(object = seurat_obj)

# 3. Generate and Save Plot
message(glue("Generating QC Violin Plots..."))
p1 <- VlnPlot(
  object   = seurat_obj,
  features = c("nCount_peaks", "nFeature_peaks", "TSS.enrichment", "nucleosome_signal"),
  pt.size  = 0.1,
  ncol     = 4
)

plot_path <- path(plot_dir, "03_atac_qc_violins.pdf")
ggsave(filename = plot_path, plot = p1, width = 14, height = 6)
message(glue("-> Saved QC plots to: {plot_path}"))

# 4. Save QC Object
save_path <- path(obj_dir, "03_atac_qc_calculated.rds")
saveRDS(seurat_obj, file = save_path)
message(glue("-> Saved QC-calculated object to: {save_path}"))
message(glue("==============================================="))

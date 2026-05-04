#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fs)
library(glue)
library(readr)
library(dplyr)
library(ggplot2)
library(future)

# Parrallel Processing (Increase threads(workers) to 8)
plan(multicore, workers = 8)
options(future.globals.maxSize = 8 * 1024^3)

message("=========Parrallel processing enabled=========")

# Fetch configurations
obj_dir    <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir   <- Sys.getenv("PLOT_DIR", unset = "results/plots")
tables_dir <- Sys.getenv("TABLES_DIR", unset = "results/tables")

message("===============================================")
message(" Step 09: ChromVAR Motif Activity Analysis     ")
message("===============================================")

# 1. Load Object
# (Loading the object from Step 08 that already contains the motif matrix)
message("Loading motif-annotated ATAC object...")
atac_obj <- readRDS(path(obj_dir, "08_atac_motifs.rds"))

# 2. Run ChromVAR
message("Calculating single-cell motif activity via ChromVAR...")
message("(This computes GC-bias backgrounds across all cells)")
atac_obj <- RunChromVAR(
  object = atac_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# 3. Switch active assay to ChromVAR
DefaultAssay(atac_obj) <- 'chromvar'

# 4. Plot Single-Cell Motif Activity
# Target Motif: MA0050.2 (JASPAR ID for IRF1 - Monocyte master regulator)
target_motif <- "MA0050.2"

message("Plotting single-cell motif activity on UMAP...")
p1 <- FeaturePlot(
  object = atac_obj,
  features = target_motif,
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.5
) + 
  ggtitle(glue("ChromVAR Activity: {target_motif} (IRF1)")) + 
  theme_minimal()

plot_path <- path(plot_dir, glue("09_chromVAR_{target_motif}_umap.pdf"))
ggsave(filename = plot_path, plot = p1, width = 8, height = 6)
message(glue("-> Saved ChromVAR UMAP plot to: {plot_path}"))

# 5. Save Final Object
saveRDS(atac_obj, path(obj_dir, "09_atac_chromvar.rds"))
message("===============================================")
message("Step 09: ChromVAR Activity Completed")

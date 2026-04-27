#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fs)
library(glue)
library(readr)
library(dplyr)
library(ggplot2)
library(future)

# Parrallel Processing (Increase threads(workers) to 8)
plan("multicore", workers = 8)
options(future.globals.maxSize = 8000 * 1024^2)

message("=========Parrallel processing enabled=========")

# Fetch configurations
obj_dir    <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir   <- Sys.getenv("PLOT_DIR", unset = "results/plots")
tables_dir <- Sys.getenv("TABLES_DIR", unset = "results/tables")

message("===============================================")
message(" Step 08: Motif Enrichment Analysis            ")
message("===============================================")

# 1. Load Object & Marker Peaks
atac_obj <- readRDS(path(obj_dir, "06_atac_annotated.rds"))
da_peaks <- read_csv(path(tables_dir, "07_marker_peaks.csv"), show_col_types = FALSE)

# 2. Extract JASPAR Motifs
message("Extracting known transcription factor motifs from JASPAR2020...")
pfm <- GetMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# 3. Add Motifs to ATAC Object
# This step maps the A/C/T/G sequences of our peaks using the hg38 genome
message("Scanning peaks for motif sequences (using hg38 genome)...")
atac_obj <- AddMotifs(
  object = atac_obj,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# 4. Find Enriched Motifs for a Specific Lineage
# Let's isolate the biologically significant peaks for CD14+ Monocytes
mono_peaks <- da_peaks %>%
  filter(cluster == "CD14+ Monocytes" & p_val_adj < 0.05) %>%
  pull(gene) # Signac stores peak coordinates in the 'gene' column

message("Testing for enriched motifs in CD14+ Monocytes...")
enriched_motifs <- FindMotifs(
  object = atac_obj,
  features = mono_peaks
)

# Save the enriched motifs table
write_csv(enriched_motifs, path(tables_dir, "08_CD14_Mono_enriched_motifs.csv"))

# 5. Plot the Top 4 Motifs
message("Plotting top 4 enriched motifs...")
top_motifs <- head(rownames(enriched_motifs), 4)
p1 <- MotifPlot(
  object = atac_obj,
  motifs = top_motifs,
  assay = 'peaks'
) + theme_minimal()

plot_path <- path(plot_dir, "08_CD14_Mono_Motifs.pdf")
ggsave(filename = plot_path, plot = p1, width = 8, height = 4)
message(glue("-> Saved Motif plot to: {plot_path}"))

# 6. Save Final Object
saveRDS(atac_obj, path(obj_dir, "08_atac_motifs.rds"))
message("===============================================")

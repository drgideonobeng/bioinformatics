#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(fs)
  library(glue)
})

message("--- Step 2: Model Execution ---")

out_dir <- Sys.getenv("OUT_DIR")
if(out_dir == "") stop("Error: OUT_DIR is empty. Sourcing config.sh failed.")

rds_in <- path(out_dir, "dds_unprocessed.rds")
if (!file_exists(rds_in)) stop(glue("Error: Cannot find {rds_in}"))

dds <- readRDS(rds_in)

message("Running DESeq2 algorithm...")
dds <- DESeq(dds)

message("Calculating Variance Stabilizing Transformation (VST) for QC...")
vsd <- vst(dds, blind = FALSE)

# Checkpoint Saves
saveRDS(dds, path(out_dir, "dds_analyzed.rds"))
saveRDS(vsd, path(out_dir, "vsd_transformed.rds"))
message(glue("=> Saved analyzed dds and vsd objects to {out_dir}/"))


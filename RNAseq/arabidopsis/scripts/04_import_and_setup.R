#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 1: Import & Setup ---")

# Paths and Metadata
quant_dir <- path("results", "quants")
out_dir <- path("results")
dir_create(out_dir)

sample_ids <- c("DRR016125", "DRR016126", "DRR016127", "DRR016128")
files <- path(quant_dir, paste0(sample_ids, "_quant"), "quant.sf")
names(files) <- sample_ids

if(!all(file_exists(files))) stop("Error: Not all quant.sf files were found.")

coldata <- data.frame(
  row.names = sample_ids,
  condition = factor(c("Control", "Control", "Treated", "Treated"))
)

# Transcript-to-Gene Mapping
message("Building tx2gene mapping...")
tmp_quant <- read_tsv(files[1], show_col_types = FALSE)
tx2gene <- data.frame(
  TXNAME = tmp_quant$Name, 
  GENEID = sub("\\..*", "", tmp_quant$Name)
)

# Import and Initialize
message("Running tximport...")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

message("Initializing DESeq2 Dataset...")
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Checkpoint Save
rds_out <- path(out_dir, "dds_unprocessed.rds")
saveRDS(dds, rds_out)
message(glue("=> Saved setup object to: {rds_out}"))

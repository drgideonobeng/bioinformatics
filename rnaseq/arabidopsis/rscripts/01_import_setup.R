#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 1: Import & Setup (Fully Dynamic) ---")

# 1. Read dynamic paths from config.sh
out_dir <- Sys.getenv("OUT_DIR")
if(out_dir == "") stop("Error: OUT_DIR is empty. Make sure run_r_pipeline.sh sources config.sh!")
dir_create(out_dir)

metadata_file <- Sys.getenv("METADATA_FILE")
if(!file_exists(metadata_file)) stop(glue("Error: Cannot find metadata file at {metadata_file}"))

# 2. Load Experimental Design
message("Loading experimental design from metadata.csv...")
coldata <- read_csv(metadata_file, show_col_types = FALSE)

# Convert to standard dataframe and set row names
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$Sample

# Dynamically set Control level
control_level <- Sys.getenv("CONTROL_TREATMENT")
if(control_level == "") stop("Error: CONTROL_TREATMENT is empty.")

coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Treatment <- relevel(coldata$Treatment, ref = control_level)

message("Experimental Design Loaded:")
print(coldata)

# 3. Build Paths to quant.sf files dynamically
sample_ids <- rownames(coldata)
# Point directly to the quants folder inside your dynamic project output directory!
quant_dir <- path(Sys.getenv("OUT_DIR"), "quants")
files <- path(quant_dir, sample_ids, "quant.sf")

if(!all(file_exists(files))) {
  missing <- files[!file_exists(files)]
  stop(glue("Error: Missing quant.sf files for: {paste(names(missing), collapse=', ')}"))
}

# 4. Transcript-to-Gene Mapping
message("\nBuilding tx2gene mapping...")
tmp_quant <- read_tsv(files[1], show_col_types = FALSE)
tx2gene <- data.frame(
  TXNAME = tmp_quant$Name, 
  GENEID = sub("\\..*", "", tmp_quant$Name)
)

# 5. Import and Initialize DESeq2
message("Running tximport...")
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

message("\nConfirming rownames of coldata match colnames of count matrix...")
if (!all(rownames(coldata) == colnames(txi$counts))) {
  stop("Error: rownames of coldata and colnames of txi$counts do not match")
} else {
  message("Success: rownames and colnames are perfectly aligned.")
}

message("\nInitializing DESeq2 Dataset...")
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ Treatment)

# Pre-filtering: remove genes with very low read counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6. Checkpoint Save
rds_out <- path(out_dir, "dds_unprocessed.rds")
saveRDS(dds, rds_out)
message(glue("=> Saved setup object to: {rds_out}"))


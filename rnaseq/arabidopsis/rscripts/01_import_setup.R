#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 1: Import & Setup (Dynamic Metadata) ---")

# Define and create the output directory so we can save at the end!
out_dir <- path("results")
dir_create(out_dir)

# 1. Read the metadata.csv file dynamically
metadata_file <- path("metadata.csv")
if(!file_exists(metadata_file)) stop("Error: Cannot find metadata.csv in the root folder!")

message("Loading experimental design from metadata.csv...")
coldata <- read_csv(metadata_file, show_col_types = FALSE)

# Convert to a standard dataframe and set the Sample column as rownames (required by DESeq2)
coldata <- as.data.frame(coldata)
rownames(coldata) <- coldata$Sample

# Ensure Treatment is a factor, and set the Control as the reference level
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Treatment <- relevel(coldata$Treatment, ref = "Water")

message("Experimental Design Loaded:")
print(coldata)

# 2. Build Paths to quant.sf files using the Sample IDs from the metadata
sample_ids <- rownames(coldata)
quant_dir <- path("results", "quants")
files <- path(quant_dir, sample_ids, "quant.sf")
names(files) <- sample_ids

if(!all(file_exists(files))) stop("Error: Not all quant.sf files were found. Did Salmon finish?")

# 3. Transcript-to-Gene Mapping
message("\nBuilding tx2gene mapping...")
tmp_quant <- read_tsv(files[1], show_col_types = FALSE)
tx2gene <- data.frame(
  TXNAME = tmp_quant$Name, 
  GENEID = sub("\\..*", "", tmp_quant$Name)
)

# 4. Import and Initialize
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

# 5. Checkpoint Save
rds_out <- path(out_dir, "dds_unprocessed.rds")
saveRDS(dds, rds_out)
message(glue("=> Saved setup object to: {rds_out}"))

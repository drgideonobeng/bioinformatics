#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(fs)
  library(glue)
})

message("--- Step 1: Import & Setup (featureCounts Pipeline) ---")

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

# 3. Read the featureCounts matrix
message("\nLoading featureCounts matrix...")
# Pointing to the counts directory defined in your config
counts_dir <- Sys.getenv("COUNTS_DIR")
count_file <- path(counts_dir, "featureCounts_matrix.txt")

if(!file_exists(count_file)) stop(glue("Error: Cannot find counts at {count_file}"))

# Read the file. featureCounts has 1 header line we skip with comment.char = "#"
# We set row.names = 1 to make the Gene IDs (Column 1) the row names
counts_raw <- read.table(count_file, header = TRUE, row.names = 1, comment.char = "#", stringsAsFactors = FALSE)

# featureCounts outputs Chr, Start, End, Strand, Length in columns 1-5. 
# We only want the actual sample count data (column 6 to the end)
count_matrix <- counts_raw[, 6:ncol(counts_raw)]

# 4. Clean up the Column Names
# R will ingest the long file paths as column names (e.g., X.Users.go140...SRR32395006_Aligned...).
# This loop extracts just the SRR ID to match the metadata.
message("Cleaning up column names to match metadata...")
colnames(count_matrix) <- sapply(colnames(count_matrix), function(col_name) {
  match <- rownames(coldata)[sapply(rownames(coldata), function(id) grepl(id, col_name))]
  if(length(match) > 0) return(match[1]) else return(col_name)
})

# 5. Import and Initialize DESeq2
message("\nConfirming rownames of coldata match colnames of count matrix...")
# Reorder count_matrix columns to strictly match the row order of coldata
count_matrix <- count_matrix[, rownames(coldata)]

if (!all(rownames(coldata) == colnames(count_matrix))) {
  stop("Error: rownames of coldata and colnames of count_matrix do not match.")
} else {
  message("Success: rownames and colnames are perfectly aligned.")
}

message("\nInitializing DESeq2 Dataset...")
dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                              colData = coldata, 
                              design = ~ Treatment)

# Pre-filtering: remove genes with fewer than 10 total reads across all samples combined
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 6. Checkpoint Save
rds_out <- path(out_dir, "dds_unprocessed.rds")
saveRDS(dds, rds_out)
message(glue("=> Saved setup object to: {rds_out}"))

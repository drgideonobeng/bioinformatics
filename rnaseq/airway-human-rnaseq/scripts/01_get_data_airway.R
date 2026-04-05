#!/usr/bin/env Rscript

# --- Step 0: Set a CRAN Mirror (Fixes the Mirror Error) ---
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- Step 1: Load/Install Libraries ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# BiocManager handles Bioconductor packages differently than CRAN
if (!require("airway", quietly = TRUE)) BiocManager::install("airway", update = FALSE)

suppressPackageStartupMessages({
library(airway)
library(tidyverse)
library(fs)
})

cat("=== Step 1: Extracting Data(Airway) ===\n\n")

# Load the built-in dataset
data(airway)

# --- Step 2: Export the Count Matrix ---
# assay(airway) gets the matrix. We turn it into a dataframe and 
# move the Ensembl Gene IDs from rownames to a real column.
counts <- as.data.frame(assay(airway)) %>%
  rownames_to_column("gene_id")

# --- Safety Check: Create folders if they don't exist ---
dir_create("data/raw")
dir_create("metadata")

write_csv(counts, "data/raw/airway_counts.csv")

# --- Step 3: Export the Metadata ---
# We use the 'Run' column for sample_id and 'dex' for the condition.
metadata <- as.data.frame(colData(airway)) %>%
  select(sample_id = Run, condition = dex) %>%
  mutate(condition = factor(condition, levels = c("untrt", "trt"), labels = c("Control", "Treated")))

write_csv(metadata, "metadata/airway_metadata.csv")

# --- Step 4: Final Verification ---
# A pro move: checking that colnames of counts (minus gene_id) match metadata rows
count_cols <- colnames(counts)[-1] # Remove 'gene_id'
metadata_ids <- metadata$sample_id

if (all(count_cols == metadata_ids)) {
  cat("[OK] Verification Success: Counts and Metadata are perfectly aligned!\n")
} else {
  cat("[WARN] Verification Failed: Sample order does not match.\n")
}

cat("[OK] Created airway_counts.csv and airway_metadata.csv.\n")

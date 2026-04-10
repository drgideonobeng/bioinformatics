#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

cat("=== Step 2: DESeq2 Analysis (TCGA-COAD) ===\n")

# 1. Load Data with check.names=FALSE to preserve the UUID hyphens
counts <- read.csv("data/cancer_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/cancer_metadata.csv", row.names = 1, check.names = FALSE)

# 2. Handle Duplicate Gene IDs
cat("-> Collapsing duplicate gene IDs...\n")
counts$base_id <- gsub("\\..*$", "", rownames(counts))
counts_collapsed <- counts %>%
  group_by(base_id) %>%
  summarise(across(everything(), sum), .groups = 'drop') %>%
  column_to_rownames("base_id")

# 3. Precision ID Syncing
cat("-> Syncing Metadata and Count IDs...\n")
# Ensure both use hyphens (UUIDs)
colnames(counts_collapsed) <- gsub("\\.", "-", colnames(counts_collapsed))
rownames(metadata) <- gsub("\\.", "-", rownames(metadata))

common_samples <- intersect(colnames(counts_collapsed), rownames(metadata))
counts_collapsed <- counts_collapsed[, common_samples]
metadata <- metadata[common_samples, ]

# 4. Labeling using the column found in your snapshot
cat("-> Labeling Tumor vs Normal...\n")
# The snapshot showed this exact name:
target_col <- "tcga.cgc_sample_sample_type"

if (target_col %in% colnames(metadata)) {
    metadata$condition <- ifelse(grepl("Tumor", metadata[[target_col]], ignore.case = TRUE), "Tumor", 
                          ifelse(grepl("Normal", metadata[[target_col]], ignore.case = TRUE), "Normal", NA))
} else {
    # Fallback to barcode logic if the column is missing for some samples
    cat("   Note: Target column missing, using barcode decoder fallback.\n")
    sample_codes <- substr(metadata$tcga.tcga_barcode, 14, 15)
    metadata$condition <- ifelse(as.numeric(sample_codes) < 10, "Tumor", "Normal")
}

# Remove any samples that couldn't be classified
metadata <- metadata[!is.na(metadata$condition), ]
counts_collapsed <- counts_collapsed[, rownames(metadata)]

metadata$condition <- factor(metadata$condition, levels = c("Normal", "Tumor"))

# Print a summary so you know it found both groups
stats <- table(metadata$condition)
cat(paste0("   Summary: Found ", stats["Normal"], " Normal and ", stats["Tumor"], " Tumor samples.\n"))

if(stats["Normal"] == 0) stop("CRITICAL ERROR: No 'Normal' samples found. Check your metadata labeling!")

# 5. Pre-filtering
cat("-> Filtering low-count genes...\n")
keep <- rowSums(counts_collapsed >= 10) >= 10
counts_filtered <- counts_collapsed[keep, ]

# 6. Run DESeq2
cat("-> Starting DESeq2 (Processing 546 samples - this will take a few minutes)...\n")
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts_filtered)),
                              colData = metadata,
                              design = ~ condition)

saveRDS(dds, "results/dds_output.rds")

dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "Tumor", "Normal"))

# 7. Save Results
if(!dir.exists("results")) dir.create("results")
write.csv(as.data.frame(res), "results/colon_cancer_de_results.csv")
cat("[OK] Analysis complete! Results saved to results/colon_cancer_de_results.csv\n")

#!/usr/bin/env Rscript

# --- Step 3.1: Load the necessary tools ---
# We use 'suppressPackageStartupMessages' to keep our logs clean.
suppressPackageStartupMessages({
  library(DESeq2)    # The engine for differential expression
  library(tidyverse) # For data manipulation (reading files, pipes)
})

cat("=== RNA-seq Step 3: Differential Expression Analysis ===\n\n")

# --- Step 3.2: Load and Prepare the Data ---
# DESeq2 requires specific formats. 
# 1. The Count Matrix needs Gene IDs as row names.
# 2. The Metadata needs Sample IDs as row names.

cat("[INFO] Loading data files...\n")
counts <- read_csv("data/raw/practice_counts_matrix.csv", show_col_types = FALSE) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

metadata <- read_csv("metadata/practice_sample_metadata.csv", show_col_types = FALSE) %>%
  column_to_rownames("sample_id")

# --- Step 3.3: Set the Comparison Baseline ---
# We must tell R which group is the 'Control'. 
# This ensures that 'Positive Log Fold Change' means 'Higher in Treatment'.
metadata$condition <- factor(metadata$condition, levels = c("Control", "Treatment"))

# --- Step 3.4: Create the DESeq2 Object ---
# This object (dds) is like a container that holds your counts, metadata, and the 
# mathematical formula (~ condition) we want to test.
cat("[INFO] Creating DESeq2 object...\n")
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

# --- Step 3.5: Run the Analysis Manually (The "Small Data" Fix) ---
cat("[INFO] Running DESeq2 steps manually for practice data...\n")

# 1. Estimate sequencing depth (size factors)
dds <- estimateSizeFactors(dds)

# 2. Estimate gene-wise noise (dispersions)
dds <- estimateDispersionsGeneEst(dds)

# 3. FORCE the final noise estimates to be the gene-wise ones
# This is the specific fix the error requested
dispersions(dds) <- mcols(dds)$dispGeneEst

# 4. Run the statistical test (Wald Test)
dds <- nbinomWaldTest(dds)

cat("[OK] Model fitting complete.\n")
# --- Step 3.6: Get Results with 'apeglm' Shrinkage ---
# In pharma, we use 'shrinkage' to stabilize the results of genes with low counts.
# This prevents 'false positives' that look like big changes but are just noise.
cat("[INFO] Extracting and shrinking results...\n")
res <- lfcShrink(dds, coef = "condition_Treatment_vs_Control", type = "apeglm")

# --- Step 3.7: Save the Results ---
# We convert the results back to a standard table (data frame) and save to CSV.
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id")

write_csv(res_df, "results/de_results.csv")
cat("[OK] Analysis complete! Results saved to 'results/de_results.csv'.\n")

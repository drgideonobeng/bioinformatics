#!/usr/bin/env Rscript

# scripts/03_annotate_results.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

cat("=== Step 3: Annotating Ensembl IDs with Gene Symbols ===\n")

# 1. Load results
if(!file.exists("results/colon_cancer_de_results.csv")) {
    stop("Error: Results file not found. Run Step 2 first!")
}
res_df <- read.csv("results/colon_cancer_de_results.csv", row.names = 1)

# 2. Clean Ensembl IDs (remove version numbers like .1, .2)
cat("-> Mapping IDs...\n")
ensembl_ids <- gsub("\\..*$", "", rownames(res_df))

# 3. Use Bioconductor to map symbols and gene names
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = ensembl_ids,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

res_df$entrez <- mapIds(org.Hs.eg.db,
                        keys = ensembl_ids,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# 4. Save the annotated results
write.csv(as.data.frame(res_df), "results/colon_cancer_annotated.csv")
cat("[OK] Annotated results saved to results/colon_cancer_annotated.csv\n")

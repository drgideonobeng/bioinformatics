#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

cat("=== Step 13: Identifying PC1 Driver Genes ===\n")

# 1. Load the VST-transformed data
if (!file.exists("results/dds_output.rds")) {
  stop("Missing dds_output.rds! Please run Step 2 first.")
}
dds <- readRDS("results/dds_output.rds")
vsd <- vst(dds, blind = FALSE)

# 2. Perform PCA manually
rv <- rowVars(assay(vsd))
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca <- prcomp(t(assay(vsd)[select_genes, ]))

# 3. Extract Loadings for PC1 and convert rownames to a column
loadings <- as.data.frame(pca$rotation) %>%
  select(PC1) %>%
  rownames_to_column("ensembl_id")

# 4. Load Annotated Results carefully
# We use row.names = 1 to handle the IDs, then move them to a column
res_annotated <- read.csv("results/colon_cancer_annotated.csv", row.names = 1) %>%
  rownames_to_column("ensembl_id") %>%
  select(ensembl_id, symbol)

# 5. Join and Find Top Drivers
pc1_top_genes <- loadings %>%
  inner_join(res_annotated, by = "ensembl_id") %>%
  mutate(absolute_loading = abs(PC1)) %>%
  arrange(desc(absolute_loading)) %>%
  head(20)

# 6. Print and Save
cat("\nTop 10 Genes Driving 24% of Sample Variance (PC1):\n")
print(head(pc1_top_genes %>% select(symbol, PC1), 10))

write.csv(pc1_top_genes, "results/pc1_driver_genes.csv", row.names = FALSE)
cat("\n[OK] Driver genes saved to results/pc1_driver_genes.csv\n")

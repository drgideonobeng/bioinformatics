#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(DESeq2)
  library(fs)
})

cat("=== Step 7: Generating Expression Heatmap ===\n")

# 1. Load the DESeq2 object
if (!file_exists("results/dds_airway.rds")) {
    stop("Error: results/dds_airway.rds not found. Please run Step 2 first!")
}
dds <- readRDS("results/dds_airway.rds")

# 2. Transform the data (VST is essential for heatmaps)
vsd <- vst(dds, blind = FALSE)

# 3. Load results and get Top 50 genes
res <- read_csv("results/airway_de_results_annotated.csv", show_col_types = FALSE)

top_genes_df <- res %>%
  filter(padj < 0.05, !is.na(symbol)) %>%
  slice_min(order_by = padj, n = 50)

# 4. Create the matrix
mat <- assay(vsd)[top_genes_df$gene_id, ]
rownames(mat) <- top_genes_df$symbol

# 5. Smart Metadata Selection
# This grabs the actual metadata table from your dds object
metadata <- as.data.frame(colData(dds))

# We'll look for 'dex' first, but fallback to whatever the 1st column is if 'dex' is missing
target_col <- if("dex" %in% colnames(metadata)) "dex" else colnames(metadata)[1]
df_ann <- metadata[, target_col, drop = FALSE]

cat(paste0("[INFO] Using '", target_col, "' for heatmap annotation.\n"))

# 6. Plotting with Z-scores
dir_create("figures")
pheatmap(mat, 
         annotation_col = df_ann, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_colnames = FALSE,
         scale = "row",         # Crucial: calculates (Value - Mean) / SD
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Top 50 DE Genes (Z-score Normalized)",
         filename = "figures/airway_heatmap.png",
         width = 8, 
         height = 10)

cat("[OK] Heatmap saved to figures/airway_heatmap.png\n")


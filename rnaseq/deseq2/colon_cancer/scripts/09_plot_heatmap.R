#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
})

cat("=== Step 9: Generating Heatmap ===\n")

# 1. Load data
counts <- read.csv("data/cancer_counts.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("data/cancer_metadata.csv", row.names = 1, check.names = FALSE)
res <- read.csv("results/colon_cancer_annotated.csv", row.names = 1)

# 2. Sync IDs and Get Top 50 Genes
top_genes <- res %>% arrange(padj) %>% head(50) %>% pull(symbol)
# Filter counts and match column names to metadata rows
colnames(counts) <- gsub("\\.", "-", colnames(counts))
rownames(metadata) <- gsub("\\.", "-", rownames(metadata))
common <- intersect(colnames(counts), rownames(metadata))

# 3. Prepare Plotting Matrix (Z-score normalization)
# We log-transform and then Z-score to show relative changes
plot_mat <- counts[rownames(res)[1:50], common] %>% as.matrix()
plot_mat <- log2(plot_mat + 1)
plot_mat <- t(scale(t(plot_mat))) # Z-score by row
rownames(plot_mat) <- res$symbol[1:50]

# 4. Define Annotation (The Color Bar at the top)
# We need to pull the 'condition' logic from Step 2
sample_codes <- substr(metadata[common, "tcga.tcga_barcode"], 14, 15)
cond <- ifelse(as.numeric(sample_codes) < 10, "Tumor", "Normal")
anno <- HeatmapAnnotation(Status = cond, 
                          col = list(Status = c("Tumor" = "firebrick", "Normal" = "dodgerblue")))

# 5. Save Heatmap
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
pdf("plots/top50_heatmap.pdf", width = 12, height = 10)
Heatmap(plot_mat, name = "Z-score", 
        top_annotation = anno,
        col = col_fun,
        show_column_names = FALSE, 
        cluster_columns = TRUE,
        column_title = "TCGA-COAD: Top 50 DE Genes")
dev.off()

cat("[OK] Heatmap saved to plots/top50_heatmap.pdf\n")

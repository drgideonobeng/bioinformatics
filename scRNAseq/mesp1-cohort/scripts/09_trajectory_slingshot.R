#!/usr/bin/env Rscript

# ========== 0. LOAD PACKAGES ==========
suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(viridis)
  library(fs)
  library(glue)
  library(readr)
  library(dplyr)
  library(future)
})

# ========== 1. SETUP & PARALLELIZATION ==========
obj_dir   <- Sys.getenv("OBJ_DIR", unset = "results/objects")
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

n_cores <- as.numeric(Sys.getenv("THREADS", unset = 8))
plan("multicore", workers = n_cores)

# ========== 2. LOAD DATA ==========
load_path <- path(obj_dir, "05_seurat_annotated.rds")
message(glue("Loading Seurat object from: {load_path}"))
seurat_obj <- readRDS(load_path)
DefaultAssay(seurat_obj) <- "RNA"

# ========== 3. SUBSET & DYNAMIC OUTLIER REMOVAL ==========
message("Subsetting Endothelial cells...")

# Tell Seurat to use your 'CellType' column for identities
Idents(seurat_obj) <- "CellType"

# Now it will successfully find and subset the Endothelial cells
cluster_03 <- subset(seurat_obj, idents = "Endothelial")

message("Removing spatial outliers...")
umap_coords <- Embeddings(cluster_03, reduction = "umap")

# Using -7 as a threshold to keep the main bottom cluster 
y_threshold <- -7 
cells_to_keep <- rownames(umap_coords)[umap_coords[, 2] < y_threshold]

if (length(cells_to_keep) == 0) {
  message(glue("Warning: No cells found with UMAP_2 < {y_threshold}."))
  message("Automatic Fallback: Keeping all Endothelial cells to prevent crash.")
  cells_to_keep <- rownames(umap_coords)
}

cluster_03 <- subset(cluster_03, cells = cells_to_keep)
cluster_03 <- JoinLayers(cluster_03) # Seurat v5 fix

# ========== 4. RUN SLINGSHOT ==========
message("Converting to SingleCellExperiment and running Slingshot...")
sce <- as.SingleCellExperiment(cluster_03, assay = "RNA")

# Running trajectory math
sce <- slingshot(sce, 
                 clusterLabels = 'seurat_clusters', 
                 reducedDim = 'UMAP')

# Add pseudotime back to Seurat metadata
cluster_03$pseudotime <- sce$slingPseudotime_1

# ========== 5. FIND CORRELATED GENES ==========
message("Calculating Spearman correlation for gene discovery...")
expr_matrix <- GetAssayData(cluster_03, assay = "RNA", layer = "data")
pt_values   <- cluster_03$pseudotime

# Only test genes expressed in at least 10% of these cells to save time
num_cells <- ncol(cluster_03)
genes_to_test <- rownames(cluster_03)[Matrix::rowSums(expr_matrix > 0) > (num_cells * 0.1)]
expr_dense <- as.matrix(expr_matrix[genes_to_test, ])

cor_results <- cor(t(expr_dense), pt_values, method = "spearman", use = "pairwise.complete.obs")

cor_df <- data.frame(
  Gene = rownames(cor_results), 
  Spearman_Rho = cor_results[, 1]
) %>%
  filter(!is.na(Spearman_Rho)) %>%
  arrange(desc(abs(Spearman_Rho)))

dir_create(table_dir)
write_csv(cor_df, path(table_dir, "07_pseudotime_correlated_genes.csv"))

# ========== 6. VISUALIZATION ==========
dir_create(plot_dir)
message("Generating trajectory plots...")

# Plot 1: Feature Plot
pdf(path(plot_dir, "09_slingshot_pseudotime.pdf"), width = 7, height = 6)
p1 <- FeaturePlot(cluster_03, features = "pseudotime") + 
      scale_color_viridis_c(option = "magma") +
      ggtitle("Endothelial Maturation: E8.0 to E10.5 Pseudotime")
print(p1)
invisible(dev.off())

# Plot 2: Trendlines for Top 6 Genes
top_genes <- head(cor_df$Gene, 6)
pdf(path(plot_dir, "09_top_cluster03_trajectory_trendlines.pdf"), width = 8, height = 6)

for (gene in top_genes) {
  plot_data <- FetchData(cluster_03, vars = c("pseudotime", gene))
  colnames(plot_data) <- c("Pseudotime", "Expression")
  
  p <- ggplot(plot_data, aes(x = Pseudotime, y = Expression)) +
       geom_point(alpha = 0.2, color = "grey70", size = 0.8) +
       geom_smooth(method = "loess", color = "firebrick", fill = "lightpink") +
       theme_minimal() +
       labs(title = glue("{gene} expression over Maturation"),
            subtitle = glue("Spearman Rho: {round(cor_df$Spearman_Rho[cor_df$Gene == gene], 3)}"),
            x = "Pseudotime", y = "Normalized Expression")
  print(p)
}
invisible(dev.off())

# ========== 7. SAVE OBJECT ==========
message("Saving final trajectory object...")
saveRDS(cluster_03, path(obj_dir, "09_Trajectory.rds"))
message("Step 09 Complete!")

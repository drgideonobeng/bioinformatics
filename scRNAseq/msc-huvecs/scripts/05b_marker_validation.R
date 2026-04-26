#!/usr/bin/env Rscript

# scripts/05b_marker_validation.R

library(Seurat)
library(tidyverse)
library(fs)
library(glue)

# 1. Environment variables
target_cluster <- Sys.getenv("CLUSTER_NAME")
obj_dir        <- Sys.getenv("OBJ_DIR",   unset = "results/objects")
plot_dir       <- Sys.getenv("PLOT_DIR",  unset = "results/plots/validation")
table_dir      <- Sys.getenv("TABLE_DIR", unset = "results/tables")

# 2. Extract Validation Genes 
# Read the file, filter to the target row, and extract the string directly
genes_string <- read_csv("validation_genes.csv", show_col_types = FALSE) %>%
  filter(cluster == target_cluster) %>%
  pull(genes)

if (length(genes_string) == 0) {
  stop(glue("Cluster '{target_cluster}' not found in validation_genes.csv"))
}

# Split the string and clean up spaces
validation_genes <- genes_string %>% 
  str_split_1(";") %>% 
  str_trim()

# 3. Load and Join Layers
message(glue("Processing {target_cluster}..."))
seurat_obj <- readRDS(path(obj_dir, "05_seurat_annotated.rds")) %>% 
  JoinLayers()

# 4. Statistical Validation (Method 1)
message("Calculating differential expression...")
stats <- FindMarkers(seurat_obj, ident.1 = target_cluster, 
                     features = validation_genes, only.pos = TRUE)
write_csv(stats, path(table_dir, glue("05b_{target_cluster}_stats.csv")))

# 5. Visualizations
message("Generating plots...")

# A. DotPlot (Summary)
p_dot <- DotPlot(seurat_obj, features = validation_genes) + RotatedAxis()
ggsave(path(plot_dir, glue("05b_DotPlot_{target_cluster}.pdf")), p_dot, width = 10, height = 6)

# B. Violin Plot (Distribution)
p_vln <- VlnPlot(seurat_obj, features = validation_genes, pt.size = 0)
ggsave(path(plot_dir, glue("05b_Violin_{target_cluster}.pdf")), p_vln, width = 12, height = 8)

# C. FeaturePlots (Merged UMAP for all samples on a single page)
# Filter out any missing genes first so the plot doesn't crash
valid_features <- validation_genes[validation_genes %in% rownames(seurat_obj)]

if (length(valid_features) > 0) {
  # Passing all genes directly creates a single grid containing all valid markers
  p_feat <- FeaturePlot(seurat_obj, features = valid_features, ncol = 2) & 
            theme(legend.position = "right")
  
  # Save the grid as a single page PDF
  ggsave(path(plot_dir, glue("05b_FeaturePlots_{target_cluster}.pdf")), p_feat, width = 12, height = 10)
  
} else {
  message("Warning: None of the validation genes were found in the Seurat object. Skipping FeaturePlot.")
}

message(glue("Successfully validated {target_cluster}!"))
message("\n=========Step 05b: Marker Validation Complete========")

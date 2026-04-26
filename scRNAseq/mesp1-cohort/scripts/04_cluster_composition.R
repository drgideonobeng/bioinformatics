#!/usr/bin/env Rscript

# scripts/04_cluster_composition.R

library(Seurat)
library(tidyverse)
library(fs)
library(glue)

# 1. Setup
obj_dir   <- Sys.getenv("OBJ_DIR")
plot_dir  <- Sys.getenv("PLOT_DIR")
seurat_obj <- readRDS(path(obj_dir, "02_seurat_harmony_integrated.rds"))

# 2. Composition Barplot (Sample vs Cluster)
message("Generating cluster composition plot...")
comp_data <- seurat_obj@meta.data %>%
  group_by(seurat_clusters, orig.ident) %>%
  tally() %>%
  group_by(seurat_clusters) %>%
  mutate(pct = n/sum(n) * 100)

p1 <- ggplot(comp_data, aes(x = seurat_clusters, y = pct, fill = orig.ident)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cluster Composition by Sample", y = "Percentage of Cells", x = "Cluster")

# 3. FeaturePlot of the "Transition" genes
message("Plotting marker expression...")
p2 <- FeaturePlot(seurat_obj, 
                  features = c("Hand1", "Nppa", "Myl2", "Nkx2-6"), 
                  ncol = 2, reduction = "umap") # Adjust reduction name if needed

# 4. Save
pdf(path(plot_dir, "04_cluster_analysis.pdf"), width = 10, height = 8)
print(p1)
print(p2)
dev.off()

message("Plots saved to results/plots/04_cluster_analysis.pdf")
message("\n=========Step 04: Cluster Analysis Complete========")

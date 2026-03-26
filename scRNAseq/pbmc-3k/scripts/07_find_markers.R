# scripts/07_find_markers.R
library(Seurat)
library(fs)
library(dplyr)
library(ggplot2)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
root_dir <- Sys.getenv("ROOT_DIR")

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")
if (root_dir == "") root_dir <- getwd()

# Ensure we have a place to save data tables
table_dir <- path(root_dir, "results", "tables")
dir_create(table_dir)

# 2. Load the clustered object
load_path <- path(obj_dir, "03_pbmc_clustered.rds")
message("Loading clustered data from: ", load_path)
pbmc <- readRDS(load_path)

# 3. Find Marker Genes
# only.pos = TRUE means we only care about genes that are turned UP in a cluster
message("Finding marker genes for all clusters... (This may take 1-2 minutes)")
pbmc.markers <- FindAllMarkers(pbmc, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)

# 4. Preview the Top 2 markers per cluster in the terminal
message("-------------------------------------------")
message("Top 2 Marker Genes per Cluster:")
top2 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
print(top2[, c("cluster", "gene", "avg_log2FC", "p_val_adj")])
message("-------------------------------------------")

# 5. Save the complete list to a CSV for biologists to review
csv_path <- path(table_dir, "01_cluster_markers.csv")
write.csv(pbmc.markers, file = csv_path, row.names = FALSE)
message("Saved full marker table to: ", csv_path)

# 6. Generate Feature Plots for canonical immune markers
# MS4A1 = B-cells, LYZ = Monocytes, CD3E = T-cells, NKG7 = NK cells
message("Generating FeaturePlots for canonical genes...")
p1 <- FeaturePlot(pbmc, 
                  features = c("MS4A1", "LYZ", "CD3E", "NKG7"), 
                  ncol = 2)

pdf_path <- path(plot_dir, "03_marker_feature_plots.pdf")
pdf(pdf_path, width = 10, height = 8)
print(p1)
dev.off()

message("Success! Marker feature plots saved to: ", pdf_path)

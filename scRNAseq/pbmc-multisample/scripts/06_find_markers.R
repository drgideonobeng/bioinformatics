
# scripts/06_find_markers.R
library(Seurat)
library(dplyr)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 06...")
obj_dir   <- Sys.getenv("OBJ_DIR")
table_dir <- Sys.getenv("TABLE_DIR")
top_genes <- as.numeric(Sys.getenv("TOP_GENES"))

# ========== 2. LOAD CLUSTERED OBJECT ==========
message("Loading clustered UMAP dataset...")
umap_obj <- readRDS(path(obj_dir, "05_umap.rds"))

# ========== 3. JOIN LAYERS (Seurat v5 Fix) ==========
message("Joining Seurat v5 data layers...")
# This glues the isolated Untreated and Treated counts together so FindAllMarkers can do math on them
umap_obj <- JoinLayers(umap_obj)

# ========== 4. FIND CLUSTER MARKERS ==========
message("Running FindAllMarkers to identify cell types...")

markers <- FindAllMarkers(umap_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# ========== 5. ORGANIZE AND SAVE RESULTS ==========
message("Organizing top markers per cluster...")
# Group the results by cluster and pull the top 5 most defining genes for each cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = top_genes, order_by = avg_log2FC)

message("Saving full marker list to results/tables...")
write.csv(markers, path(table_dir, "06_cluster_markers.csv"), row.names = FALSE)

message("Saving top 5 marker list to results/tables...")
write.csv(top_markers, path(table_dir, "06a_cluster_top_markers.csv"), row.names = FALSE)

message("Showing top 5 most defining genes for each cluster")
print(top_markers, n = Inf)
message("Step 06 Complete! Biomarkers identified and saved.")


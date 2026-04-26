# scripts/08_find_markers.R
library(Seurat)
library(fs)
library(dplyr)
library(glue)

obj_dir <- Sys.getenv("OBJ_DIR")
if (obj_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 1. Load the clustered object
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 2. Find all markers
message("Calculating marker genes for all clusters...")
message("(This requires a lot of math and may take a couple of minutes!)")

# only.pos = TRUE means we only care about genes that are turned UP in a cluster
markers <- FindAllMarkers(seurat_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

# 3. Save the results to a CSV file
tables_dir <- path("results", "tables")
dir_create(tables_dir)
save_path <- path(tables_dir, "08_cluster_markers.csv")

write.csv(markers, file = save_path, row.names = FALSE)

# 4. Print a sneak peek to the terminal
message("\n--- TOP 5 MARKER GENES PER CLUSTER ---")
top_genes <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n = 5, order_by = avg_log2FC)

# 5. Save the sheet of top genes
top_genes_path <- path(tables_dir, "08_top_cluster_markers_cheatsheet.csv")
write.csv(top_genes, file = top_genes_path, row.names = FALSE)
 
print(top_genes[, c("cluster", "gene", "avg_log2FC")], n= Inf)

message("\nStep 08 Complete! Full marker list saved to: ", save_path)
message("Top genes cheat sheet saved to: ", top_genes_path)

# scripts/08a_annotate_manual.R
library(Seurat)
library(ggplot2)
library(fs)

message("Loading configuration variables for Step 08a (Manual)...")
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

message("Loading clustered dataset...")
umap_obj <- readRDS(path(obj_dir, "05_umap.rds"))

message("Manually assigning biological names to clusters...")
# Our custom dictionary!
cluster_names <- c(
  "0" = "CD14+ Monocytes",
  "1" = "CD4 T-cells",
  "2" = "CD4 T-cells", 
  "3" = "CD8 T-cells",
  "4" = "CD16+ Monocytes",
  "5" = "B-cells (Stimulated)",
  "6" = "B-cells (Resting)",
  "7" = "Dendritic Cells",
  "8" = "NK Cells",
  "9" = "Plasmacytoid DCs",
  "10" = "Megakaryocytes",
  "11" = "Activated T-cells",
  "12" = "Unknown 12",
  "13" = "Unknown 13",
  "14" = "Unknown 14",
  "15" = "Red Blood Cells"
)

# Apply the manual names
umap_obj <- RenameIdents(umap_obj, cluster_names)
umap_obj$manual_cell_type <- Idents(umap_obj)

message("Generating manual annotated UMAP...")
pdf(path(plot_dir, "08a_manual_UMAP.pdf"), width = 10, height = 6)
print(DimPlot(umap_obj, reduction = "umap", group.by = "manual_cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Manual Annotation: Human Detective Work"))
dev.off()

message("Saving manual annotated Seurat object...")
saveRDS(umap_obj, path(obj_dir, "08a_annotated.rds"))
message("Step 08a Complete!")

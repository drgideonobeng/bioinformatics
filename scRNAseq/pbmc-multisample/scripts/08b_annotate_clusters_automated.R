# scripts/08b_annotate_automated.R
library(Seurat)
library(ggplot2)
library(fs)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

message("Loading configuration variables for Step 08b (Automated)...")
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

message("Loading clustered dataset...")
# FIXED: Using fs::path to prevent Bioconductor masking
umap_obj <- readRDS(fs::path(obj_dir, "05_umap.rds"))
umap_obj <- JoinLayers(umap_obj) # Safety first for v5!

message("Downloading Human Immune Reference Database...")
ref_data <- celldex::MonacoImmuneData()

message("Converting Seurat object for SingleR...")
sce_obj <- as.SingleCellExperiment(umap_obj)

message("Running SingleR automated annotation (This may take a minute!)...")
predictions <- SingleR(test = sce_obj, assay.type.test = 1, ref = ref_data, labels = ref_data$label.main)

message("Transferring AI labels back to Seurat object...")
umap_obj$automated_cell_type <- predictions$labels
Idents(umap_obj) <- "automated_cell_type"

message("Generating automated annotated UMAP...")
# FIXED: Using fs::path
pdf(fs::path(plot_dir, "08b_automated_UMAP.pdf"), width = 10, height = 6)
print(DimPlot(umap_obj, reduction = "umap", group.by = "automated_cell_type", label = TRUE, repel = TRUE) +
  ggtitle("Automated Annotation: SingleR + Monaco Reference"))
dev.off()

message("Saving automated annotated Seurat object...")
# FIXED: Using fs::path
saveRDS(umap_obj, fs::path(obj_dir, "08b_annotated.rds"))
message("Step 08b Complete!")

#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(ggplot2)
library(fs)
library(glue)

# Fetch configurations
obj_dir <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")
matrix_dir <- Sys.getenv("DATA_MATRIX_DIR", unset = "data/raw_matrix")
dim_start <- as.numeric(Sys.getenv("LSI_DIMS_START", unset = 2))
dim_end <- as.numeric(Sys.getenv("LSI_DIMS_END", unset = 30))

message(glue("==============================================="))
message(glue(" Step 06: scATAC-seq Label Transfer            "))
message(glue("==============================================="))

# 1. Load the clustered ATAC object
read_path <- path(obj_dir, "05_atac_clustered.rds")
message(glue("Loading clustered ATAC object: {read_path}"))
atac_obj <- readRDS(read_path)

# 2. Build the "Translator" (Gene Activity Matrix)
message(glue("Calculating Gene Activity Matrix (This acts as our bridge to RNA)..."))
gene_activities <- GeneActivity(atac_obj)

# Add the gene activity matrix to the Seurat object as a new assay
atac_obj[['RNA']] <- CreateAssayObject(counts = gene_activities)
atac_obj <- NormalizeData(
  object = atac_obj,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac_obj$nCount_RNA)
)

# 3. Download & Load the scRNA-seq Reference Data
ref_path <- path(matrix_dir, "pbmc_rna_reference.rds")
if (!file_exists(ref_path)) {
  message(glue("Downloading annotated scRNA-seq PBMC reference..."))
  download.file(
    url = "https://www.dropbox.com/s/3f3p5nxrn5b3y4y/pbmc_10k_v3.rds?dl=1",
    destfile = ref_path,
    mode = "wb",
    quiet = TRUE
  )
}

message(glue("Loading RNA reference object..."))
rna_obj <- readRDS(ref_path)
rna_obj <- UpdateSeuratObject(rna_obj) # Ensure compatibility with Seurat v5

# 4. Find Anchors between ATAC and RNA
message(glue("Finding mathematical anchors between ATAC and RNA datasets..."))
transfer_anchors <- FindTransferAnchors(
  reference = rna_obj,
  query = atac_obj,
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

# 5. Transfer the Labels!
message(glue("Transferring Cell Type Labels..."))
predicted_labels <- TransferData(
  anchorset = transfer_anchors,
  refdata = rna_obj$celltype, # This is the metadata column in the RNA object containing the names
  weight.reduction = atac_obj[['lsi']],
  dims = dim_start:dim_end
)

# Add the predicted labels to our ATAC metadata
atac_obj <- AddMetaData(object = atac_obj, metadata = predicted_labels)

# 6. Plot the Final Annotated UMAP
message(glue("Plotting Annotated UMAP..."))
p1 <- DimPlot(atac_obj, group.by = "predicted.id", label = TRUE, pt.size = 0.5) +
  ggtitle("scATAC-seq PBMCs (Annotated via Label Transfer)") +
  theme_minimal() +
  NoLegend()

plot_path <- path(plot_dir, "06_atac_annotated_umap.pdf")
ggsave(filename = plot_path, plot = p1, width = 8, height = 6)
message(glue("-> Saved Annotated UMAP plot to: {plot_path}"))

# 7. Save Final Object
save_path <- path(obj_dir, "06_atac_annotated.rds")
saveRDS(atac_obj, file = save_path)
message(glue("-> Saved fully annotated object to: {save_path}"))
message(glue("==============================================="))

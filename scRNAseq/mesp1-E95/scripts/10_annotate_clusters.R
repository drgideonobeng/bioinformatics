# scripts/10_annotate_clusters.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories and variables from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
project_name <- Sys.getenv("PROJECT_NAME") # Fix for the missing plot title!

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "04_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Define the biological identities based on our marker genes
message("Translating cluster numbers to biological cell types...")
cluster_identities <- c(
  "0" = "CD4+ T Cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD16+ Monocytes",
  "3" = "T Cells (FHIT+)",
  "4" = "Naive B Cells",
  "5" = "MAIT Cells",
  "6" = "NK Cells",
  "7" = "Plasma Cells",
  "8" = "Cytotoxic CD8+ T Cells",
  "9" = "Dendritic Cells (LYPD2+)",
  "10" = "CD8+ T Cells",
  "11" = "T/NK Subset",
  "12" = "Platelets (GP9+)",
  "13" = "AXL+ Dendritic Cells",
  "14" = "Memory B Cells",
  "15" = "pDCs (LRRC26+)",
  "16" = "B Cells (JCHAIN+)",
  "17" = "Basophils (CPA3+)",
  "18" = "Unknown/Rare (KIF3C+)"
)

# 4. Apply the new names to the object
seurat_obj <- RenameIdents(seurat_obj, cluster_identities)
seurat_obj$cell_type <- Idents(seurat_obj)

# 5. Generate the final Annotated UMAP Plot
message("Generating final annotated UMAP plot...")
annotated_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle(paste(project_name, "Annotated Cell Types")) +
  theme_bw() +
  theme(legend.position = "right")

# Save the plot
dir_create(plot_dir)
pdf_path <- path(plot_dir, "05_umap_annotated.pdf")
pdf(pdf_path, width = 10, height = 7) 
print(annotated_umap)
dev.off()

# 6. Save the final object
final_obj_path <- path(obj_dir, "05_seurat_annotated.rds")
saveRDS(seurat_obj, file = final_obj_path)

message("Pipeline Complete! Annotated UMAP saved to: ", pdf_path)

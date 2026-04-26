#!/usr/bin/env Rscript

# scripts/05_cluster_annotation.R

library(ggplot2)
library(Seurat)
library(fs)

# 1. Setup paths
obj_dir <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load the integrated object
message("Loading integrated Seurat object...")
seurat_obj <- readRDS(path(obj_dir, "02_seurat_harmony_integrated.rds"))

# 3. Create your Annotation Dictionary
message("Annotating clusters...")
cluster_ids <- c(
  "0"  = "Pharyngeal_Mesoderm",   # Osr1+, Foxf1+
  "1"  = "Mesoderm_Progenitors",   # Pax9+, Foxd1+
  "2"  = "Skeletal_Muscle",        # Myf5+ (Skeletal lineage)
  "3"  = "Endothelial",            # Cldn5+, Cdh5+ (Blood vessels)
  "4"  = "Somitic_Mesoderm",       # Meox1+, Pax3+
  "5"  = "Epithelium",             # Krt8+
  "6"  = "Lateral_Plate_Mesoderm", # Gsc+, Alx3+
  "7"  = "Erythroid",              # Hbq1b+, Gypa+ (Red blood cells)
  "8"  = "Progenitors_Unspec",
  "9"  = "Early_Cardiomyocytes",   # Transition state (Nkx2-6+)
  "10" = "Posterior_Mesoderm",     # Hox10/11 genes
  "11" = "Cardiomyocytes",         # YOUR TARGET (Nppa+, Myl2+)
  "12" = "Cardiac_Fibroblasts",
  "13" = "Epicardium",             # Aldh1a1+, Upk1b+
  "14" = "Myeloid_Immune"          # C1qa+, Pf4+ (Macrophages)
)

# Apply the names to the object
seurat_obj <- RenameIdents(seurat_obj, cluster_ids)
# Save it as a formal metadata column
seurat_obj$CellType <- Idents(seurat_obj)

# 4. Generate an Annotated UMAP (merged)
message("Generating Annotated UMAP...")
pdf(path(plot_dir, "05_annotated_umap_merged.pdf"), width = 10, height = 7)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) + 
  ggtitle("Annotated Clusters (E8.0 - E10.5)")
dev.off()

# 5. Save the  annotated object
message("Saving annotated object...")
save_path <- path(obj_dir, "05_seurat_annotated.rds")
saveRDS(seurat_obj, save_path)
message("Success! Object saved to: ", save_path)

message("\n=========Step 05: Cluster Annotation Complete========")

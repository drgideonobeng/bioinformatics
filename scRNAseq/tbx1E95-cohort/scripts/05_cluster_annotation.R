#!/usr/bin/env Rscript

# scripts/05_cluster_annotation.R

library(ggplot2)
library(Seurat)
library(fs)
library(glue)

# 1. Setup paths
obj_dir <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")
umap_title <- Sys.getenv("UMAP_TITLE", unset = "E9.5")

# 2. Load the integrated object
message("Loading integrated Seurat object...")
seurat_obj <- readRDS(path(obj_dir, "02_seurat_integrated.rds"))

# 3. Create your Annotation Dictionary
message("Annotating clusters...")
cluster_ids <- c(
  "0"  = "Paraxial_Mesoderm",             # Meox1+, Cyp26a1+, Sfrp1+
  "1"  = "Mesenchyme_Fibroblasts",        # Foxf2+, Lum+, Vtn+, Cthrc1+
  "2"  = "Cycling_Progenitors",           # Nusap1+, Prc1+, Kif18a+
  "3"  = "Cranial_Neural_Crest",          # Msx1+, Msx2+, Zic1+, Dkk1+
  "4"  = "Cranial_Mesoderm",              # Pitx2+, Msc+, Lhx2+
  "5"  = "Pharyngeal_Endoderm",           # Pax1+, Foxi2+, Krt7+
  "6"  = "Posterior_Pharyngeal_Endoderm", # Hoxb1+, Hoxb3+, Pax8+
  "7"  = "Splanchnic_Mesoderm",           # Hand2+, Gata5+, Alx1+
  "8"  = "Endothelial",                   # Cdh5+, Cldn5+, Icam2+, Esam+
  "9"  = "Posterior_SHF",                 # Tbx5+, Aldh1a2+, Foxf1+, Osr1+
  "10" = "Cranial_Ganglia_Neuroblasts",   # Neurod1+, Phox2b+, Tlx3+
  "11" = "Cardiomyocytes",                # Myh6+, Myl2+, Actn2+, Myl3+
  "12" = "Neural_Tube_CNS",               # Pax6+, Pax3+, Hes5+, Rfx4+
  "13" = "Cranial_Placodes_Otic"          # Sox10+, Hmx2+, Hmx3+, Dlx5+
)

# Apply the names to the object
seurat_obj <- RenameIdents(seurat_obj, cluster_ids)
# Save it as a formal metadata column
seurat_obj$CellType <- Idents(seurat_obj)

# 4. Generate an Annotated UMAP (merged)
message("Generating Annotated UMAP...")
pdf(path(plot_dir, "05_annotated_umap_merged.pdf"), width = 10, height = 7)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) + 
  ggtitle(glue("{umap_title} Annotated Clusters"))
invisible(dev.off())

# 5. Save the  annotated object
message("Saving annotated object...")
save_path <- path(obj_dir, "05_seurat_annotated.rds")
saveRDS(seurat_obj, save_path)
message("Success! Object saved to: ", save_path)

message("\n=========Step 05: Cluster Annotation Complete========")

# scripts/10_annotate_clusters.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories and variables from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")
project_name <- Sys.getenv("PROJECT_NAME")

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Define the biological identities based on our marker genes
message("Translating cluster numbers to biological cell types...")
cluster_identities <- c(
  "0" = "Posterior Arch Mesenchyme (Hox+)", 
  "1" = "Second Heart Field (Aldh1a2+)", 
  "2" = "Otic Ectoderm (Pax8+)", 
  "3" = "Cardiac Pacemaker (Shox2+)", 
  "4" = "Cranial Mesoderm (Tlx1+)", 
  "5" = "Proximal Arch / Cardiac Prog. (Isl1+)", 
  "6" = "Pharyngeal Ectoderm (Nkx2-3+)", 
  "7" = "Erythrocytes (Hba-x+)", 
  "8" = "Endothelial (Cldn5+)", 
  "9" = "Epithelium (Aldh1a1+)", 
  "10" = "Cardiomyocytes (Myh7+)", 
  "11" = "Platelets (Pf4+)"
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

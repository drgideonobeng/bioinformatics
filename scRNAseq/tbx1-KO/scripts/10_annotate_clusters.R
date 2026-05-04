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
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Define the biological identities based on our KO E9.5 marker genes
message("Translating cluster numbers to KO biological cell types...")
cluster_identities <- c(
  "0" = "Pharyngeal Mesenchyme",
  "1" = "Cranial Neural Crest",
  "2" = "Posterior Endoderm / Mesoderm",
  "3" = "Hindbrain",
  "4" = "Proliferating / Cycling Cells",
  "5" = "Transitional Mesenchyme",
  "6" = "Pharyngeal Mesoderm",
  "7" = "Anterior SHF / Arch Mesenchyme",
  "8" = "Pharyngeal Epithelium",
  "9" = "Posterior SHF / Epicardial Progenitors",
  "10" = "Endothelium",
  "11" = "Cranial Ganglia / Neuroblasts",
  "12" = "Cardiomyocytes"
)

# 4a. Apply the new names to the active identities
seurat_obj <- RenameIdents(seurat_obj, cluster_identities) 

# 4b. Save these annotations into the metadata for easy access later
seurat_obj$cell_type <- Idents(seurat_obj)

# 5. Generate the Annotated UMAP Plot
message("Generating annotated UMAP plot...")
annotated_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle(paste(project_name, "Annotated Cell Lineages")) +
  theme_bw() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))

# 6. Save the plot
dir_create(plot_dir)
pdf_path <- path(plot_dir, "10_umap_annotated.pdf")
pdf(pdf_path, width = 10, height = 7) 
print(annotated_umap)
invisible(dev.off())


# 7. Save the final object
message("Saving annotated Seurat object...")
final_obj_path <- path(obj_dir, "10_seurat_annotated.rds")
saveRDS(seurat_obj, file = final_obj_path)

message("Step 10 Complete! Annotated UMAP saved to: ", pdf_path)

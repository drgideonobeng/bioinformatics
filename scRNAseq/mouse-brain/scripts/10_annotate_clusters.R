# scripts/10_annotate_clusters.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "04_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Define the biological identities based on our marker genes
message("Translating cluster numbers to biological cell types...")
cluster_identities <- c(
  "0" = "Nos1+ Neurons",
  "1" = "Layer 6 Projection Neurons",
  "2" = "Dividing Neural Progenitors",
  "3" = "SST+ Interneurons",
  "4" = "Intermediate Progenitors (IPCs)",
  "5" = "Slit3+ Specialized Cells",
  "6" = "HTR3A+ Interneurons",
  "7" = "Layer 5 Projection Neurons"
)

# 4. Apply the new names to the object
seurat_obj <- RenameIdents(seurat_obj, cluster_identities)

# We also save this newly named identity into a specific metadata column for safekeeping
seurat_obj$cell_type <- Idents(seurat_obj)

# 5. Generate the final Annotated UMAP Plot
message("Generating final annotated UMAP plot...")
annotated_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle("Developing Mouse Cortex: Annotated Cell Types") +
  theme_bw() +
  theme(legend.position = "right")

# Save the plot
dir_create(plot_dir)
pdf_path <- path(plot_dir, "05_umap_annotated.pdf")
pdf(pdf_path, width = 10, height = 7) # Made slightly wider to fit the legend!
print(annotated_umap)
dev.off()

# 6. Save the final, fully processed and annotated Seurat object
final_obj_path <- path(obj_dir, "05_seurat_annotated.rds")
saveRDS(seurat_obj, file = final_obj_path)

message("Pipeline Complete! Annotated UMAP saved to: ", pdf_path)
message("Final Seurat object saved to: ", final_obj_path)

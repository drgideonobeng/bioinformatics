# scripts/08_annotate_cells.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Setup paths
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "03_pbmc_clustered.rds")
message("Loading clustered data...")
pbmc <- readRDS(load_path)

# 3. Define the cell type identities
# This maps the classic PBMC 3k clusters to their biological names
new.cluster.ids <- c(
  "0" = "Naive CD4 T", 
  "1" = "Memory CD4 T", 
  "2" = "CD14+ Mono", 
  "3" = "B", 
  "4" = "CD8 T", 
  "5" = "FCGR3A+ Mono", 
  "6" = "NK", 
  "7" = "DC", 
  "8" = "Platelet"
)

# 4. Apply the new names to the object
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# 5. Generate the Final Annotated UMAP
message("Generating annotated UMAP...")
p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() +
  ggtitle("Annotated PBMC Cell Types") +
  theme_bw()

# Save the plot
dir_create(plot_dir)
pdf_path <- path(plot_dir, "04_final_annotated_umap.pdf")
pdf(pdf_path, width = 8, height = 6)
print(p1)
dev.off()

# 6. Save the final, presentation-ready Seurat Object
save_path <- path(obj_dir, "04_pbmc_annotated.rds")
saveRDS(pbmc, file = save_path)

message("--------------------------------------------------")
message("Pipeline Complete!")
message("Final Annotated Object: ", save_path)
message("Final UMAP Plot: ", pdf_path)
message("--------------------------------------------------")

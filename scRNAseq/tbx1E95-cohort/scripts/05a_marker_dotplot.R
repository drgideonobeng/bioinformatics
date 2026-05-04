#!/usr/bin/env Rscript

# scripts/05a_marker_dotplot.R

library(Seurat)
library(ggplot2)
library(fs)

# 1. Setup paths from your config logic
obj_dir   <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load the Annotated Object (from script 05)
# If you haven't run script 05 yet, use "02_seurat_integrated.rds"
load_path <- path(obj_dir, "05_seurat_annotated.rds")
if (!file_exists(load_path)) {
    message("Annotated object not found, falling back to integrated object...")
    load_path <- path(obj_dir, "02_seurat_integrated.rds")
}
seurat_obj <- readRDS(load_path)

# 3. Define Representative Genes for the Global View
# These are selected from your 03_top_markers.csv to create a clear "Diagonal"
genes_to_plot <- c(
  "Osr1", "Foxf1",             # Cluster 0: Pharyngeal
  "Pax9", "Foxd1",             # Cluster 1: Progenitors
  "Myf5", "Vgll3",             # Cluster 2: Skeletal Muscle
  "Cldn5", "Cdh5",             # Cluster 3: Endothelial
  "Meox1", "Pax3",             # Cluster 4: Somitic
  "Krt8", "Alcam",             # Cluster 5: Epithelium
  "Gsc", "Alx3",               # Cluster 6: Lateral Plate
  "Hbq1b", "Gypa",             # Cluster 7: Erythroid
  "Rpgrip1", "Notch2",         # Cluster 8: Unspec
  "Nkx2-6", "Mylk3",           # Cluster 9: Early Cardiomyocytes
  "Hoxc10", "Hoxa11",          # Cluster 10: Posterior
  "Nppa", "Myl2", "Tnnt2",     # Cluster 11: Cardiomyocytes (TARGET)
  "Papss2", "Erbb3",           # Cluster 12: Fibroblasts
  "Aldh1a1", "Upk1b",          # Cluster 13: Epicardium
  "C1qa", "Pf4"                # Cluster 14: Myeloid
)

# 4. Generate the DotPlot
message("Generating Global DotPlot...")
p <- DotPlot(seurat_obj, 
             features = unique(genes_to_plot), 
             cols = c("lightgrey", "blue"), # Blue is standard, or use "red" for heat
             dot.scale = 8) + 
  RotatedAxis() +
  theme(axis.text.x = element_text(face = "italic", size = 10)) +
  labs(title = "Global Marker Signature: Mesp1 Lineage Atlas",
       subtitle = "Dot size = % of cells expressing | Color = Expression Intensity",
       x = "Marker Gene",
       y = "Cluster / Cell Type")

# 5. Save the Plot
save_path <- path(plot_dir, "05a_global_marker_dotplot.pdf")
pdf(save_path, width = 14, height = 8)
print(p)
dev.off()

message("Success! Global View saved to: ", save_path)

message("\n=========Step 05a: Marker Dotplot Complete========")

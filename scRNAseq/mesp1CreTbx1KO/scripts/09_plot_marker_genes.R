# scripts/09_plot_marker_genes.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "" || plot_dir == "") stop("Environment variables missing. Run 'source config.sh' first.")

# 2. Load the clustered object
load_path <- path(obj_dir, "06_seurat_clustered.rds")
message("Loading clustered data from: ", load_path)
seurat_obj <- readRDS(load_path)

# 3. Define the biological marker genes we want to map
features_to_plot <- c("Ppbp", "Cldn5", "Cxcl16", "Nppa")

message("Generating FeaturePlots...")

# 4. Generate the plot
# 'order = TRUE' ensures the colored dots are drawn on TOP of the gray dots!
p <- FeaturePlot(seurat_obj, 
                 features = features_to_plot, 
                 pt.size = 0.5, 
                 order = TRUE, 
                 ncol = 2)

# 5. Save the plot to a PDF
dir_create(plot_dir)
save_path <- path(plot_dir, "04_marker_feature_plots.pdf")

pdf(save_path, width = 10, height = 8)
print(p)
dev.off()

message("Step 09 Complete! Feature plots saved to: ", save_path)

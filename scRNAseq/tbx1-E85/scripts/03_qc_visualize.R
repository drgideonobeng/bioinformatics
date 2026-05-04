# scripts/03_qc_visualize.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the unfiltered object safely using fs::path
message("Loading unfiltered Seurat object...")
seurat_obj <- readRDS(path(obj_dir, "02_seurat_unfiltered.rds"))

# 3. Generate Violin Plots
message("Generating QC Violin Plots...")
dir_create(plot_dir)

# 4. Save the plot
pdf(path(plot_dir, "03_qc_violins.pdf"), width = 10, height = 6)
print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
dev.off()

message("Step 03 Complete. QC plots saved to: ", plot_dir)


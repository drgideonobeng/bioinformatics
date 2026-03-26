# scripts/04_scale_and_pca.R
library(Seurat)
library(ggplot2)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 04...")
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

# ========== 2. LOAD NORMALIZED OBJECT ==========
message("Loading normalized dataset...")
norm_obj <- readRDS(path(obj_dir, "03_normalized.rds"))

# ========== 3. SCALE THE DATA ==========
message("Scaling data (shifting mean to 0, variance to 1)...")
# We apply scaling to all genes in the dataset
all_genes <- rownames(norm_obj)
scaled_obj <- ScaleData(norm_obj, features = all_genes)

# ========== 4. RUN PCA ==========
message("Running Principal Component Analysis (PCA)...")
# We run PCA only on our 2000 Highly Variable Features to save massive amounts of RAM
pca_obj <- RunPCA(scaled_obj, features = VariableFeatures(object = scaled_obj))

# ========== 5. GENERATE ELBOW PLOT ==========
message("Generating PCA Elbow Plot...")
# The Elbow plot helps us visually confirm how many PCs capture the real biology
pdf(path(plot_dir, "04_pca_elbow_plot.pdf"), width = 8, height = 6)
print(ElbowPlot(pca_obj, ndims = 50))
dev.off()

# ========== 6. SAVE PCA OBJECT ==========
message("Saving PCA-reduced Seurat object...")
saveRDS(pca_obj, path(obj_dir, "04_pca.rds"))

message("Step 04 Complete! Data scaled, PCA calculated, and Elbow plot generated.")

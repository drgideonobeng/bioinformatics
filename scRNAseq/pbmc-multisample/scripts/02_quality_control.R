# scripts/02_qc.R
library(Seurat)
library(ggplot2)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 02...")
obj_dir          <- Sys.getenv("OBJ_DIR")
plot_dir         <- Sys.getenv("PLOT_DIR")

min_genes        <- as.numeric(Sys.getenv("MIN_GENES"))
max_genes        <- as.numeric(Sys.getenv("MAX_GENES"))
max_mt           <- as.numeric(Sys.getenv("MAX_MT_PERCENT"))
mt_pattern       <- Sys.getenv("MT_PATTERN")

# ========== 2. LOAD RAW OBJECT ==========
message("Loading raw merged object...")
merged_obj <- readRDS(path(obj_dir, "01_raw_merged.rds"))

# ========== 3. QUALITY CONTROL & FILTERING ==========
message("Calculating mitochondrial percentages...")
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = mt_pattern)

message("Generating PRE-filter QC plot...")
pdf(path(plot_dir, "02a_pre_filter_QC.pdf"), width = 10, height = 5)
print(VlnPlot(merged_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "condition", pt.size = 0.1))
dev.off()

message("Applying strict config filtering thresholds...")
filtered_obj <- subset(merged_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes & percent.mt < max_mt)

message("Generating Post-filter QC plot...")
pdf(path(plot_dir, "02b_post_filter_QC.pdf"), width = 10, height = 5)
print(VlnPlot(filtered_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "condition", pt.size = 0))
dev.off()

# ========== 4. SAVE FILTERED OBJECT ==========
message("Saving filtered Seurat object...")
saveRDS(filtered_obj, path(obj_dir, "02_qc_filtered.rds"))

message("Step 02 Complete! Quality control finished.")
message("------------------------------------------------")
message(glue("Cells BEFORE filtering: {ncol(merged_obj)}"))
message(glue("Cells AFTER filtering:  {ncol(filtered_obj)}"))
message("------------------------------------------------")

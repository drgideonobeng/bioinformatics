# scripts/03_qc_visualize.R
library(Seurat)
library(ggplot2)
library(fs)

# 1. Pull directories from config
obj_dir  <- Sys.getenv("OBJ_DIR")
plot_dir <- Sys.getenv("PLOT_DIR")

if (obj_dir == "") stop("OBJ_DIR is empty. Please run 'source config.sh' first.")

# 2. Load the unfiltered object safely using fs::path
load_path <- path(obj_dir, "01_pbmc_unfiltered.rds")
pbmc <- readRDS(load_path)

# 3. Extract metadata for simple, robust plotting
md <- pbmc@meta.data

# 4. Build Independent Plots
p1 <- ggplot(md, aes(x = "PBMC", y = nFeature_RNA)) + 
  geom_violin(fill = "skyblue") + theme_bw() + ggtitle("Unique Genes")

p2 <- ggplot(md, aes(x = "PBMC", y = nCount_RNA)) + 
  geom_violin(fill = "lightgreen") + theme_bw() + ggtitle("Total Molecules")

p3 <- ggplot(md, aes(x = "PBMC", y = percent.mt)) + 
  geom_violin(fill = "salmon") + theme_bw() + ggtitle("Mitochondrial %")

s1 <- ggplot(md, aes(x = nCount_RNA, y = percent.mt)) + 
  geom_point(alpha = 0.4) + theme_bw() + ggtitle("Molecules vs. Mito %")

s2 <- ggplot(md, aes(x = nCount_RNA, y = nFeature_RNA)) + 
  geom_point(alpha = 0.4) + theme_bw() + ggtitle("Molecules vs. Unique Genes")

# 5. Save to PDF (One per page to avoid layout crashes)
dir_create(plot_dir) # Ensure the plot folder exists before saving
pdf_path <- path(plot_dir, "01_qc_pre_filter_plots.pdf")

pdf(pdf_path, width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(s1)
print(s2)
dev.off()

message("Success! QC plots saved to: ", pdf_path)

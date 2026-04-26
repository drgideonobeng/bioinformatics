#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(ggplot2)
library(fs)
library(glue)

# 1. Setup paths
obj_dir   <- Sys.getenv("OBJ_DIR", unset = "results/objects")
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load Object
message("Loading Seurat object...")
seurat_obj <- readRDS(fs::path(obj_dir, "05_seurat_annotated.rds"))

# 3. CRITICAL FOR SEURAT V5: Join layers before DE
message("Joining RNA layers...")
seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")

# 4. Isolate the Cardiomyocytes
message("Subsetting Cardiomyocytes...")
cluster10_obj <- subset(seurat_obj, idents = "Cardiomyocytes")

# 5. Set identity to the sample timepoints
# Based on your previous composition plot, the orig.ident values are Mesp1_E105 and Mesp1_E95
Idents(cluster10_obj) <- cluster10_obj$orig.ident

# 6. Run Differential Expression
message("Running DE between E105 and E95 within Cluster 10...")
# ident.1 is the "numerator" (positive Log2FC means UP in E105)
# ident.2 is the "denominator" (negative Log2FC means UP in E95)
de_results <- FindMarkers(cluster10_obj, 
                          ident.1 = "Mesp1_E105", 
                          ident.2 = "Mesp1_E95", 
                          assay = "RNA",          
                          logfc.threshold = 0.25, 
                          min.pct = 0.1)

# Format the results cleanly
de_results$gene <- rownames(de_results)
de_results <- de_results %>%
  relocate(gene) %>%
  filter(!grepl("^Hb", gene)) %>%
  arrange(p_val_adj) # Sort by significance

dir_create(table_dir)
save_path <- fs::path(table_dir, "05_Cluster10_E105_vs_E95_DE_filtered.csv")
write.csv(de_results, file = save_path, row.names = FALSE)
# ------------------------------

message(" -> Results saved to: ", save_path)
# Format the results cleanly
de_results$gene <- rownames(de_results)
de_results <- de_results %>%
  relocate(gene) %>%
  arrange(p_val_adj) # Sort by significance

message(" -> Results saved to: ", save_path)

# 7. Generate a Volcano Plot
message("Generating Volcano Plot...")
# Create a new column to categorize the genes for coloring
de_results$significance <- "Not Significant"
de_results$significance[de_results$avg_log2FC > 0.5 & de_results$p_val_adj < 0.05] <- "UP in E105"
de_results$significance[de_results$avg_log2FC < -0.5 & de_results$p_val_adj < 0.05] <- "UP in E95"

p <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("UP in E95" = "blue", "Not Significant" = "grey", "UP in E105" = "red")) +
  theme_minimal() +
  labs(title = "Cluster 10: E10.5 vs E9.5",
       subtitle = "Positive Log2FC = Enriched in E10.5 | Negative Log2FC = Enriched in E9.5",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5)

# Save Plot
dir_create(plot_dir)
pdf(path(plot_dir, "05_volcano_cluster10.pdf"), width = 8, height = 6)
print(p)
dev.off()

message("Success! Volcano plot saved to: ", path(plot_dir, "05_volcano_cluster10.pdf"))

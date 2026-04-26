library(Seurat)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(viridis)

# 1.  Set up paths
obj_dir   <- Sys.getenv("OBJ_DIR", unset = "results/objects")
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load your object with the correct path
message("Loading Seurat object...")
seurat_obj <- readRDS(fs::path(obj_dir, "05_seurat_annotated.rds"))
DefaultAssay(seurat_obj) <- "RNA" # Ensure RNA is the default assay

# 3. Subset Cardiomyocytes
message("Subsetting Cardiomyocytes...")
c10 <- subset(seurat_obj, idents = "Cardiomyocytes")

# 4. Remove outliers
message("Removing spatial outliers...")
umap_coords <- Embeddings(c10, reduction = "umap") # Extract the UMAP coordinates
# Find the names of the cells that are in the main upper cluster (umap_2 > 3)
cells_to_keep <- rownames(umap_coords)[umap_coords[, 2] > 3]
# Subset the object to only keep those good cells
c10 <- subset(c10, cells = cells_to_keep)

# 5. THE V5 FIX: Join the split layers back together
message("Joining Seurat v5 layers...")
c10 <- JoinLayers(c10)

# 6. Convert to SingleCellExperiment 
message("Converting to SingleCellExperiment format...")
sce <- as.SingleCellExperiment(c10, assay = "RNA")

# 7. Run Slingshot
message("Running Slingshot trajectory...")
sce <- slingshot(sce, 
                 clusterLabels = 'seurat_clusters', 
                 reducedDim = 'UMAP')

# 8. Extract Pseudotime
# This adds the pseudotime values to your Seurat object metadata
c10$pseudotime <- sce$slingPseudotime_1

# 9. Visualization
message("Plotting trajectory...")
plot_path <- "results/plots/07_slingshot_pseudotime.pdf"

pdf(plot_path, width = 7, height = 6)
FeaturePlot(c10, features = "pseudotime") + 
  scale_color_viridis_c(option = "magma") +
  ggtitle("Cardiomyocyte Maturation: E9.5 to E10.5 Pseudotime")
dev.off()

message("Success! Trajectory mapped and plot saved.")

# ---------------------------------------------------------
# 10. Save the Trajectory-Embedded Object
# ---------------------------------------------------------
message("Saving updated Seurat object...")
saveRDS(c10, "results/objects/07_Cluster10_Trajectory.rds")

# ---------------------------------------------------------
# 11. Find Top Pseudotime-Correlated Genes (Discovery)
# ---------------------------------------------------------
message("Running Spearman correlation across all genes...")

expr_matrix <- GetAssayData(c10, assay = "RNA", layer = "data")
pt_values <- c10$pseudotime
expr_dense <- as.matrix(expr_matrix)

cor_results <- cor(t(expr_dense), pt_values, method = "spearman", use = "pairwise.complete.obs")

cor_df <- data.frame(Gene = rownames(cor_results), Spearman_Rho = cor_results[, 1])
cor_df <- na.omit(cor_df)
cor_df <- cor_df[order(-abs(cor_df$Spearman_Rho)), ] # Sort by strongest driver

# Save the master list
write.csv(cor_df, "results/tables/07_pseudotime_correlated_genes.csv", row.names = FALSE)

# ---------------------------------------------------------
# 12. Generate Trendlines for the Top Drivers (Visualization)
# ---------------------------------------------------------
# AUTOMATION: Automatically grab the top 6 genes from the correlation results
top_dynamic_genes <- head(cor_df$Gene, 6)

message("Plotting trendlines for top dynamically expressed genes: ", paste(top_dynamic_genes, collapse = ", "))

trendline_path <- "results/plots/07_top_trajectory_trendlines.pdf"
pdf(trendline_path, width = 6, height = 5)

for (gene in top_dynamic_genes) {
  
  plot_data <- FetchData(c10, vars = c("pseudotime", gene))
  colnames(plot_data) <- c("Pseudotime", "Expression")
  
  p <- ggplot(plot_data, aes(x = Pseudotime, y = Expression)) +
    geom_point(alpha = 0.4, color = "darkgray", size = 1) +
    geom_smooth(method = "loess", color = "firebrick", fill = "lightpink") +
    theme_minimal() +
    labs(title = paste(gene, "Expression over Maturation"),
         subtitle = paste("Spearman Rho:", round(cor_df$Spearman_Rho[cor_df$Gene == gene], 3)),
         x = "Pseudotime (E9.5 -> E10.5)",
         y = "Normalized Expression")
  
  print(p)
}

dev.off()
message("Pipeline Complete! Top genes found and plotted.")

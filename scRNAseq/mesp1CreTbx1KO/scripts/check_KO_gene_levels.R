library(Seurat)
library(dplyr)
library(glue)

# 1. Load the object
seurat_obj <- readRDS("results/objects/06_seurat_clustered.rds")

# 2. Percentage Breakdown
message("\n--- Tbx1 Detection Frequency (Percentage of Cells) ---")
# a. Get a logical vector: TRUE if Tbx1 expression > 0 (FIX: Added $Tbx1)
tbx1 <- FetchData(seurat_obj, vars = "Tbx1")$Tbx1 > 0
# b. Create a contingency table vs the sample IDs
counts_table <- table(seurat_obj$orig.ident, tbx1)
# c. Calculate proportions (margin = 1 makes it row-wise percentages)
results <- prop.table(counts_table, margin = 1)
print(results)

# 3. Intensity Breakdown
message("\n--- Tbx1 Raw Count Distribution by Sample ---")
samples <- unique(seurat_obj$orig.ident)

for (s in samples) {
  message(glue("\nSample: {s}"))
  counts <- seurat_obj[["RNA"]]$counts["Tbx1", seurat_obj$orig.ident == s]
  
  if(sum(counts > 0) > 0) {
    print(table(counts[counts > 0]))
    message("Average intensity in positive cells: ", round(mean(counts[counts > 0]), 2))
  } else {
    message("No Tbx1 molecules detected.")
  }
}

# 4. Cluster Localization: Which clusters have the most Tbx1?
message("\n--- Tbx1 Distribution Across Clusters ---")
tbx1_positive_cells <- subset(seurat_obj, subset = Tbx1 > 0)
cluster_counts <- table(tbx1_positive_cells$seurat_clusters)
print(cluster_counts)

# 5. Lineage Diagnostic: Are these cells Endoderm or Ectoderm?
# We check for Sox17 (Endoderm) and Krt18 (Ectoderm)
message("\n--- Lineage Marker Co-expression Check ---")
# This creates a data frame of Tbx1+ cells and their lineage markers
lineage_check <- FetchData(tbx1_positive_cells, 
                           vars = c("Tbx1", "Sox17", "Krt18", "Mesp1", "seurat_clusters"))

# Calculate how many Tbx1+ cells are actually other lineages
endo_count <- sum(lineage_check$Sox17 > 0)
ecto_count <- sum(lineage_check$Krt18 > 0)
mesp_count <- sum(lineage_check$Mesp1 > 0)
total_tbx1 <- nrow(lineage_check)

message(glue("Total Tbx1+ cells: {total_tbx1}"))
message(glue("Tbx1+ cells co-expressing Sox17 (Endoderm): {endo_count} ({round(endo_count/total_tbx1*100, 1)}%)"))
message(glue("Tbx1+ cells co-expressing Krt18 (Ectoderm): {ecto_count} ({round(ecto_count/total_tbx1*100, 1)}%)"))
message(glue("Tbx1+ cells co-expressing Mesp1 (Mesoderm Lineage): {mesp_count} ({round(mesp_count/total_tbx1*100, 1)}%)"))

# 6. Final Visual Visuals (Saves to your plot directory)
plot_dir <- Sys.getenv("PLOT_DIR", unset = "results/plots")
pdf(file.path(plot_dir, "06_tbx1_diagnostic_plots.pdf"), width = 12, height = 10)

# FIX: Moved DotPlot here and wrapped in print()
print(DotPlot(seurat_obj, features = "Tbx1", group.by = "orig.ident") + RotatedAxis())
# Show where Tbx1 is vs where the Lineage markers are
print(FeaturePlot(seurat_obj, features = c("Tbx1", "Sox17", "Krt18", "Mesp1"), ncol = 2))
# Violin plot to see intensity per cluster
print(VlnPlot(seurat_obj, features = "Tbx1", group.by = "seurat_clusters"))

dev.off()

message("\nDiagnostic plots saved to: ", plot_dir)

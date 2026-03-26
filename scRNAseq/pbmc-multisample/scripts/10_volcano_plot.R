# scripts/10_volcano_plot.R
library(ggplot2)
library(dplyr)
library(fs)
library(ggrepel)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 10...")
table_dir <- Sys.getenv("TABLE_DIR")
plot_dir  <- Sys.getenv("PLOT_DIR")

# ========== 2. LOAD DIFFERENTIAL GENES ==========
message("Loading the 5,553 significant genes...")
# We set row.names = 1 because Seurat saves the gene names in the first unnamed column
de_results <- read.csv(fs::path(table_dir, "09_treated_vs_untreated_Bcells.csv"), row.names = 1)

# Create a proper column for gene names so ggplot can read them
de_results$gene <- rownames(de_results)

# ========== 3. SET MATH THRESHOLDS ==========
message("Categorizing genes into UP, DOWN, or Not Significant...")
# We are looking for genes that doubled in expression (log2FC > 1 or < -1) 
# AND have a significant p-value (p_val_adj < 0.05)

de_results$Significance <- "Not Significant"
de_results$Significance[de_results$avg_log2FC > 1 & de_results$p_val_adj < 0.05] <- "Upregulated (Treated)"
de_results$Significance[de_results$avg_log2FC < -1 & de_results$p_val_adj < 0.05] <- "Downregulated (Treated)"

# ========== 4. ISOLATE THE SMOKING GUNS ==========
message("Finding the top 15 most extreme genes to label...")
# We sort by the lowest p-value and highest fold-change, then grab the top 15
top_genes <- de_results %>% 
  filter(Significance != "Not Significant") %>% 
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>% 
  head(15)

# ========== 5. GENERATE VOLCANO PLOT ==========
message("Drawing Volcano Plot...")

pdf(fs::path(plot_dir, "10_Bcell_Volcano_Plot.pdf"), width = 10, height = 8)

volcano <- ggplot(de_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated (Treated)" = "blue", 
                                "Not Significant" = "grey80", 
                                "Upregulated (Treated)" = "red")) +
  # Add the dashed threshold lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  # Add the labels for the top 15 genes so they don't overlap
  geom_text_repel(data = top_genes, aes(label = gene), color = "black", box.padding = 0.5) +
  theme_minimal() +
  labs(title = "Volcano Plot: B-cells (Treated vs. Untreated)",
       subtitle = "Red = Turned ON by Treated | Blue = Turned OFF by Treated",
       x = "Log2 Fold Change (Effect Size)",
       y = "-Log10 Adjusted P-Value (Statistical Certainty)")

print(volcano)
dev.off()

message("Step 10 Complete! The Volcano has erupted.")

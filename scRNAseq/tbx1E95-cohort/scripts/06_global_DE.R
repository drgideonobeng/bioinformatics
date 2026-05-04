#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)
library(fs)
library(glue)

# 1. Set up the paths
obj_dir   <- Sys.getenv("OBJ_DIR",   unset = "results/objects")
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR",  unset = "results/plots")
cluster_dir <- dir_create(table_dir, "cluster_dir")

# 2. Load the fully integrated object
message("Loading Seurat object...")
seurat_obj <- readRDS(path(obj_dir, "05_seurat_annotated.rds"))
DefaultAssay(seurat_obj) <- "RNA"

# 3. Join layers to prevent Seurat v5 matrix errors during DE testing
message("Joining layers...")
seurat_obj <- JoinLayers(seurat_obj)

# 4. Identify all annotated clusters in the dataset
clusters <- levels(Idents(seurat_obj))

# Initialize an empty list to store summary results (much faster than base R rbind!)
summary_list <- list()

message("Beginning Global Differential Expression (E10.5 vs E8.0) across all clusters...")

# 5. Loop through every cluster
for (cluster in clusters) {
  message(glue("Testing Cluster {cluster} ..."))

  # tryCatch handles clusters missing from one of the timepoints gracefully
  de_results <- tryCatch({
    FindMarkers(seurat_obj, 
                subset.ident = cluster,      
                group.by = "orig.ident",     # Update if your metadata column is named differently
                ident.1 = "Mesp1_E105", 
                ident.2 = "Mesp1_E80", 
                min.pct = 0.25,
                logfc.threshold = 1.0)
  }, error = function(e) {
    message(glue("  -> Skipping Cluster {cluster}: Not enough cells in both timepoints."))
    return(NULL)
  })
  
  # 6. Process and save the results (if the test was successful)
  if (!is.null(de_results) && nrow(de_results) > 0) {
    
    # Modern tidyverse cleaning: move row names to a proper column
    de_clean <- de_results %>% 
      rownames_to_column(var = "gene")
    
    # Write the individual cluster's full DE list cleanly
    file_path <- path(cluster_dir, glue("06_Cluster_{cluster}_E105_vs_E80_DE.csv"))
    write_csv(de_clean, file_path)
    
    # Count significant genes using tidy filtering
    sig_count <- de_clean %>% 
      filter(p_val_adj < 0.05) %>% 
      nrow()
    
    # Store the result in our list
    summary_list[[cluster]] <- tibble(Cluster = cluster, Total_Significant_DEGs = sig_count)
    
  } else {
    # If skipped or no results, record 0 DEGs
    summary_list[[cluster]] <- tibble(Cluster = cluster, Total_Significant_DEGs = 0)
  }
}

# 7. Compile and sort the summary table
# bind_rows combines the list into one dataframe, arrange sorts it dynamically
summary_df <- bind_rows(summary_list) %>% 
  arrange(desc(Total_Significant_DEGs))

# 8. Save the master summary table
summary_path <- path(table_dir, "06_Global_DE_Summary.csv")
write_csv(summary_df, summary_path)

message("--------------------------------------------------")
message(glue("Step 06 : Global DE complete! Summary saved to: {summary_path}"))
message("Here are your top most dynamically changing clusters:")
print(head(summary_df, 10))

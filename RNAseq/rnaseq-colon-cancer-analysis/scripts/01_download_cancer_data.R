#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(recount3)
  library(tidyverse)
  library(fs)
})

cat("=== Searching for TCGA Colorectal Cancer (COAD) ===\n")

# 1. Get projects
human_projects <- available_projects()

# 2. Find TCGA COAD - This is the primary Colon Adenocarcinoma dataset
proj_info <- subset(
  human_projects, 
  project == "COAD" & project_type == "data_sources"
)

# 3. If COAD isn't there, let's find a large SRP project as a backup
if (nrow(proj_info) == 0) {
    cat("[INFO] TCGA COAD not found, searching for large SRA projects...\n")
    # This finds projects with > 50 samples, which are usually the good cancer ones
    proj_info <- human_projects %>% 
      filter(n_samples > 50, file_source == "sra") %>% 
      head(1)
}

cat("[INFO] Project Found! Selected ID:", proj_info$project, "\n")

# 4. Create RSE and Transform
cat("[INFO] Downloading and creating RSE object...\n")
rse <- create_rse(proj_info)

cat("[INFO] Transforming coverage to counts...\n")
# This scales base-pair coverage to estimated read counts
assay(rse, "counts") <- transform_counts(rse)

# 5. Extract and Save
counts_rounded <- round(assay(rse, "counts"))
metadata <- as.data.frame(colData(rse))

dir_create("data")
write.csv(counts_rounded, "data/cancer_counts.csv")
write.csv(metadata, "data/cancer_metadata.csv")

cat("[OK] Data successfully saved to data/ folder!\n")

#!/usr/bin/env Rscript
library(fs)
library(glue)
library(curl)

# 1. Pull base directories from config
data_matrix_dir <- Sys.getenv("DATA_MATRIX_DIR")
if (data_matrix_dir == "") stop("DATA_MATRIX_DIR is empty. Check config.sh")

# 2. Load the sample metadata
metadata_path <- "sample_metadata.csv"
if (!file_exists(metadata_path)) stop("Missing sample_metadata.csv in the root directory!")
meta <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Base GEO download URL structure
base_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc={geo_acc}&format=file&file={geo_file}"

message("Starting multi-sample download from GEO...")

# 3. Loop through each sample in the CSV
for (i in 1:nrow(meta)) {
  geo_acc <- meta$geo_accession[i]
  sid <- meta$sample_id[i]
  
  # The exact prefix used by the authors on GEO
  # Note: You may need to adjust this if different samples have different naming conventions!
  geo_pref <- glue("{geo_acc}_{sid}") 
  
  # Create a unique sub-folder for this sample
  sample_dir <- path(data_matrix_dir, sid)
  dir_create(sample_dir)
  
  files_to_download <- list(
    "barcodes.tsv.gz" = glue("{geo_pref}_barcodes.tsv.gz"),
    "features.tsv.gz" = glue("{geo_pref}_genes.tsv.gz"), 
    "matrix.mtx.gz"   = glue("{geo_pref}_matrix.mtx.gz")
  )
  
  message(glue("\nDownloading {sid} ({geo_acc}) into {sample_dir}..."))
  
  for (name in names(files_to_download)) {
    geo_file <- files_to_download[[name]]
    download_url <- glue(base_url)
    dest_file <- path(sample_dir, name)
    
    message(glue(" -> Fetching {name}..."))
    
    tryCatch({
      curl_download(url = download_url, destfile = dest_file, quiet = TRUE)
    }, error = function(e) {
      message(glue("    [ERROR] Failed to download {name}. Check GEO accession or file prefix."))
    })
  }
}

message("\nMulti-sample download complete! Ready for Harmony Integration.")

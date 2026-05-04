#!/usr/bin/env Rscript
library(fs)
library(glue)
library(purrr)

data_matrix_dir <- Sys.getenv("DATA_MATRIX_DIR", unset = "data/raw_matrix")
dir_create(data_matrix_dir)
options(timeout = 3600)

# Pull URLs dynamically from config
tenx_files <- list(
  "matrix.h5"            = Sys.getenv("URL_MATRIX"),
  "fragments.tsv.gz"     = Sys.getenv("URL_FRAGMENTS"),
  "fragments.tsv.gz.tbi" = Sys.getenv("URL_FRAG_INDEX")
)

if (any(tenx_files == "")) stop("Missing URLs in config.sh. Please run 'source config.sh' first.")

message(glue("======================================================"))
message(glue(" Initiating Download to: {data_matrix_dir}"))
message(glue("======================================================"))

# Tidyverse approach to looping through lists
walk2(names(tenx_files), tenx_files, \(dest_name, download_url) {
  dest_file_path <- path(data_matrix_dir, dest_name)
  message(glue(" -> Fetching {dest_name}... "))
  download.file(url = download_url, destfile = dest_file_path, mode = "wb", quiet = TRUE)
})

message(glue("Download Complete! H5 matrix and fragments are ready."))

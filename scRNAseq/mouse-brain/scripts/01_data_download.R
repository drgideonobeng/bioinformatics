# scripts/01_data_download.R
library(fs)

url <- Sys.getenv("RAW_DATA_URL")
root_dir <- Sys.getenv("ROOT_DIR")

if (url == "") stop("RAW_DATA_URL is empty. Check config.sh")

# Define paths dynamically
data_dir <- path(root_dir, "data")
dir_create(data_dir)
dest_file <- path(data_dir, "raw_data.tar.gz")

message("Downloading 10x dataset from: ", url)
download.file(url, destfile = dest_file, quiet = FALSE)

message("Extracting data into: ", data_dir)
untar(dest_file, exdir = data_dir)

message("Download and extraction complete!")

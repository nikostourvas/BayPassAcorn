#!/usr/bin/env Rscript

# Script to merge WZA output chunk files
library(data.table)

# Get input parameters from Snakemake
input_dir <- dirname(snakemake@input[[1]])
output_file <- snakemake@output[[1]]
# Use the wildcards directly instead of params
sample_name <- snakemake@wildcards[["sample"]]
env_factor <- snakemake@wildcards[["envfactor"]]
chunks <- snakemake@params[[1]]

# Get file type from input filename - detect whether it is processing BF or spearman
file_type <- if(grepl("spearman", basename(output_file))) "spearman" else "BF"

# Create log connection for output
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt")
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

# Log starting info
cat(paste("Starting merging of", chunks, "chunks for", sample_name, "and factor", env_factor, "\n"))
cat(paste("File type:", file_type, "\n"))

# Get list of all chunk files using the appropriate prefix
chunk_files <- paste0(input_dir, "/", sample_name, "_", env_factor, "_", file_type, "_WZA_output_", 1:chunks, ".csv")

# Check if all files exist - abort if any are missing
missing_files <- chunk_files[!file.exists(chunk_files)]
if (length(missing_files) > 0) {
  stop("ERROR: The following chunk files are missing: ", paste(missing_files, collapse = ", "))
}

# Read first file to get header
header <- readLines(chunk_files[1], n = 1)

# Open output file and write header
writeLines(header, output_file)

# Process each chunk file
processed_rows <- 0
for (i in 1:chunks) {
  chunk_file <- chunk_files[i]
  # Read file, skipping header
  chunk_data <- fread(chunk_file, skip = 1)
  # Append to output
  fwrite(chunk_data, output_file, append = TRUE, col.names = FALSE)
  processed_rows <- processed_rows + nrow(chunk_data)
  cat(paste("Processed chunk", i, "with", nrow(chunk_data), "rows\n"))
}

cat(paste("Merged", processed_rows, "total rows to", output_file, "\n"))

# Close log connection
sink(type = "message")
sink()
close(log_conn)
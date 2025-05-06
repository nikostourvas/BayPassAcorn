#!/usr/bin/env Rscript

# Script to merge WZA output chunk files
library(data.table)

# Get input parameters from Snakemake
input_dir <- dirname(snakemake@input[[1]])
output_file <- snakemake@output[[1]]
sample_name <- snakemake@params[[2]]
env_factor <- snakemake@params[[3]]
chunks <- snakemake@params[[1]]

# Create log connection for output
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt")
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

# Log starting info
cat(paste("Starting merging of", chunks, "chunks for", sample_name, "and factor", env_factor, "\n"))

# Get list of all chunk files
chunk_files <- paste0(input_dir, "/", sample_name, "_", env_factor, "_BF_WZA_output_", 1:chunks, ".csv")

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
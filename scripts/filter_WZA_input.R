library(data.table)

# Load the data
input_data <- fread(snakemake@input[[1]])

# Filter the data
input_data <- input_data[input_data$MAF >= snakemake@params[[1]],]

# Split data into user-defined chunks
nrows <- nrow(input_data)
chunk_size <- ceiling(nrows / snakemake@params[[2]])

# Create output directory if it doesn't exist
output_dir <- dirname(snakemake@output[[1]])
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Get base filename without extension
sample_name <- snakemake@params[[3]]
dir_path <- dirname(snakemake@output[[1]])

# Split and export data chunks
for (i in 1:snakemake@params[[2]]) {
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, nrows)
  
  # Skip if we've reached the end of the data
  if (start_idx > nrows) break
  
  # Extract chunk
  chunk <- input_data[start_idx:end_idx,]
  
  # Create output filename
  output_file <- file.path(dir_path, paste0(sample_name, "_WZA_input_filtered_chunk", i, ".csv"))
  
  # Write chunk to file
  fwrite(chunk, output_file)
  
  # Print progress
  message(paste("Exported chunk", i, "with", nrow(chunk), "rows to", output_file))
}
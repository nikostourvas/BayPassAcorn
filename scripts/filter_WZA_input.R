library(data.table)

# Load the data
input_data <- fread(snakemake@input[[1]])

# Filter the data
input_data <- input_data[input_data$MAF >= snakemake@params[[1]],]

# Export the data
fwrite(input_data, snakemake@output[[1]])
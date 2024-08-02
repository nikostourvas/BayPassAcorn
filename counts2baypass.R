library(dplyr)

# Read the original data
counts <- read.table(snakemake@input[[1]],
                     sep = "\t", header = TRUE)
original_data <- counts[ ,-c(1:4)]

# Replace NAs with 0. This is ok because BayPass identifies missing data as "0 0".
original_data[is.na(original_data)] <- 0

# List of unique pool IDs
pool_ids <- unique(sub("(.ref.cnt|.alt.cnt)$", "", colnames(original_data)))

# Initialize an empty data frame with the same number of rows as the original data
result <- data.frame(matrix(ncol = length(pool_ids), nrow = nrow(original_data)))

# Name the columns of the result data frame
colnames(result) <- pool_ids

# Loop through each pool ID and merge the .ref and .alt columns
for (pool_id in pool_ids) {
  result[[pool_id]] <- paste(original_data[[paste0(pool_id, ".ref.cnt")]], 
                             original_data[[paste0(pool_id, ".alt.cnt")]], 
                             sep = " ")
}

# Write the result to a new file, columns should be separated by spaces
write.table(result, snakemake@output[[1]], sep = " ", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


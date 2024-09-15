### Nikos Tourvas 2024-07-29
### Script used to generate the complementary input files for Baypass

library(dplyr)

# Regular expression to input the files containing the poolnames for each dataset
datasets = list.files(path = "data", pattern="*_poolnames")
datasets = gsub("_poolnames", "", datasets)

generate_complementary_Baypass_inputs = function(x) {
  # Load complementary data
  pops = read.table(paste0("data/", x, "_poolnames"), header = FALSE)
  envfactors = read.csv(snakemake@input[[1]], header = TRUE)
  poolsizes = read.csv(snakemake@input[[2]], header = TRUE)
  
  # Filter the rows of the 'poolsizes' file that correspond to the populations in the vector 'pops'. Pops are labeled as 'Plot_ID' in the poolsizes file
  poolsizes_sub = poolsizes %>% 
    filter(Plot_ID %in% pops)
  # Write the filtered poolsizes to a new file.
  # Only the 'Poolsize' column should be written.
  # The 'Poolsize' column should be printed in one line (i.e. transposed) with values separated with space
  write(poolsizes_sub$poolsize, paste0("data/", x, "_poolsizes"), 
              sep = " ", ncolumns=length(poolsizes_sub$poolsize))
  
  # Filter the rows of the 'envfactors' file that correspond to the populations in the vector 'pops'. Pops are labeled as 'Plot_ID' in the envfactors file
  # Also remove columns Plot_ID and Pair_ID
  envfactors_sub = envfactors %>% 
    filter(Plot_ID %in% pops) %>% 
    select(-c("Plot_ID", "Pair_ID", "Site_description",
              "Latitude", "Longitude", "elevation")) 

  # write the column names (envfactors) to a new file
  write.table(colnames(envfactors_sub), paste0("data/", x, "_efile_envfactor_names"), sep = " ", 
              col.names=FALSE, row.names=FALSE, quote=FALSE)

  # Transpose the envfactors_sub dataframe
  envfactors_sub = t(envfactors_sub)

  # Write the filtered climate/topographic data to a new file.
  write.table(envfactors_sub, paste0("data/", x, "_efile"), sep = " ", 
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # Create the ecotype (contrast) file
  envfactors_sub = envfactors %>% 
    filter(Plot_ID %in% pops) %>% 
    select(c("Site_description")) %>% 
    t()
  # Write the filtered climate/topographic data to a new file.
  write.table(envfactors_sub, paste0("data/", x, "_ecotype"), sep = " ", 
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}

# Apply the function to all datasets
lapply(datasets, generate_complementary_Baypass_inputs)


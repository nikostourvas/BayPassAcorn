### Nikos Tourvas 2024-07-29
### Custom script with limited reusability, used to generate the complementary input files for Baypass

library(dplyr)

# Regular expression to input the files containing the poolnames for each dataset
datasets = list.files(path = "data", pattern="*HDplot_poolnames")
datasets = gsub("_poolnames", "", datasets)
lapply(datasets, generate_complementary_Baypass_inputs)

generate_complementary_Baypass_inputs = function(x) {
  # Load complementary data
  pops = read.table(paste0("data/", x, "_poolnames"), header = FALSE)
  ecotype = read.csv("data/20240729_ACORN_dem_TOPO_by_POP_nikos.csv", header = TRUE)
  poolsizes = read.csv("data/20240307_poolsizes_per_pop.csv", header = TRUE)
  
  # Filter the rows of the 'poolsizes' file that correspond to the populations in the vector 'pops'. Pops are labeled as 'Plot_ID' in the poolsizes file
  poolsizes_sub = poolsizes %>% 
    filter(Plot_ID %in% pops)
  
  # Write the filtered poolsizes to a new file.
  # Only the 'Poolsize' column should be written.
  # The 'Poolsize' column should be printed in one line (i.e. transposed) with values separated with space
  write(poolsizes_sub$poolsize, paste0("data/", x, "_poolsizes"), 
              sep = " ", ncolumns=length(poolsizes_sub$poolsize))
  
  # Filter the rows of the 'ecotype' file that correspond to the populations in the vector 'pops'. Pops are labeled as 'Plot_ID' in the ecotype file
  # Also remove columns Plot_ID and Pair_ID
  ecotype_sub = ecotype %>% 
    filter(Plot_ID %in% pops) %>% 
    select(-c("Plot_ID", "Pair_ID")) %>% 
    t()

  # Write the filtered climate/topographic data to a new file.
  write.table(ecotype_sub, paste0("data/", x, "_ecotype"), sep = " ", 
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}

source("baypass_utils.R")
# all.res.core = concatenate_res(
#                         anaprefix = snakemake@params[[1]],
#                         nsubsets=20,
#                         snpdet_prefix = snakemake@params[[4]], 
#                         retrieve_pi_xtx=TRUE,
#                         retrieve_bfis=FALSE,
#                         retrieve_c2=FALSE)

# write.table(all.res.core, snakemake@output[[1]], sep = ",", 
#             row.names = FALSE, quote = FALSE)

# all.res.std = concatenate_res(
#                         anaprefix = snakemake@params[[2]],
#                         nsubsets=20,
#                         snpdet_prefix = snakemake@params[[4]], 
#                         retrieve_pi_xtx=TRUE,
#                         retrieve_bfis=TRUE,
#                         retrieve_c2=FALSE)

# write.table(all.res.std, snakemake@output[[2]], sep = ",",
#             row.names = FALSE, quote = FALSE)

all.res.contrast = concatenate_res(
                        anaprefix = snakemake@params[[1]],
                        nsubsets=snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]], 
                        retrieve_pi_xtx=TRUE,
                        retrieve_bfis=TRUE,
                        retrieve_c2=TRUE)

write.table(all.res.contrast, snakemake@output[[1]], sep = ",", 
            row.names = FALSE, quote = FALSE)
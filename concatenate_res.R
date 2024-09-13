source("baypass_utils.R")

all.res = concatenate_res(
                        anaprefix = snakemake@params[[1]],
                        nsubsets=snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]], 
                        retrieve_pi_xtx=TRUE,
                        retrieve_bfis=TRUE,
                        retrieve_c2=snakemake@params[[4]])

write.table(all.res, snakemake@output[[2]], sep = ",", 
            row.names = FALSE, quote = FALSE)
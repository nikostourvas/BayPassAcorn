# load libraries
library(poolfstat)

# generate a pooldata object
psizes <- as.numeric(read.table(snakemake@input[[2]], header = FALSE, sep = " "))
pnames <- as.character(read.table(snakemake@input[[3]], header = FALSE, sep = " "))

poolreadcount <- vcf2pooldata(vcf.file = snakemake@input[[1]], 
                              poolsizes = psizes, poolnames = pnames)

poolreadcount

# convert the pooldata object to BayPass input files (i.e. genobaypass, snpdet, poolsize)
# and create multiple subsets that cover the whole data set
pooldata2genobaypass(poolreadcount, writing.dir = "results/subsets", 
                     prefix = snakemake@params[[1]],
                     subsamplesize = nrow(poolreadcount@readcoverage)/snakemake@params[[2]],
                     subsamplingmethod = "thinning")

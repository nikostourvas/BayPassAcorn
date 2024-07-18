# load libraries
library(poolfstat)

# generate a pooldata object
#psizes <- as.numeric(c('40', '34', '40', '40', '40', '38', '40', '36', '40', '40', '40', '40', '34', '32', '40', '40', '40', '40', '40', '30', '40', '40', '40', '34', '40', '40', '40', '40', '34', '34', '40', '38', '40', '40', '34', '32', '38', '40', '40', '40'))
#pnames <- as.character(c('300', '308', '320', '364', '390', '312', '301', '309', '321', '365', '324', '302', '310', '340', '384', '366', '325', '303', '311', '341', '385', '367', '304', '342', '386', '322', '305', '313', '343', '387', '323', '306', '314', '362', '388', '360', '307', '315', '363', '389'))
psizes <- as.numeric(read.table(snakemake@params[[1]], header = FALSE, sep = " "))
pnames <- as.character(read.table(snakemake@params[[2]], header = FALSE, sep = " "))

poolreadcount <- vcf2pooldata(vcf.file = snakemake@input[[1]], 
                              poolsizes = psizes, poolnames = pnames)

poolreadcount

# convert the pooldata object to BayPass input files (i.e. genobaypass, snpdet, poolsize)
# and create multiple subsets that cover the whole data set
pooldata2genobaypass(poolreadcount, writing.dir = "results", 
                     prefix = snakemake@params[[3]],
                     subsamplesize = nrow(poolreadcount@readcoverage)/snakemake@params[[4]],
                     subsamplingmethod = "thinning")

# load libraries
library(poolfstat)
# library(xlsx)

# set working directory
setwd("/local_disks/disk_3/BayPassAcorn")

# generate a pooldata object
psizes <- as.numeric(c('40', '34', '40', '40', '40', '38', '40', '36', '40', '40', '40', '40', '34', '32', '40', '40', '40', '40', '40', '30', '40', '40', '40', '34', '40', '40', '40', '40', '34', '34', '40', '38', '40', '40', '34', '32', '38', '40', '40', '40'))
pnames <- as.character(c('300', '308', '320', '364', '390', '312', '301', '309', '321', '365', '324', '302', '310', '340', '384', '366', '325', '303', '311', '341', '385', '367', '304', '342', '386', '322', '305', '313', '343', '387', '323', '306', '314', '362', '388', '360', '307', '315', '363', '389'))

poolreadcount <- vcf2pooldata(vcf.file = "Data/ACORN_VCF_Qrobur_minCov20IndelFilteredSNPs_Biallelic_miss0.1_ADP104_Subsampled_0.01.vcf", 
                              poolsizes = psizes, poolnames = pnames)
poolreadcount
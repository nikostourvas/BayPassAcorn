# source baypass utils functions
source("baypass_utils.R")

#Compare the estimates of Omega obtained with different subsets of SNPs
omega1 <- as.matrix(read.table(snakemake@input[[1]]))
omega2 <- as.matrix(read.table(snakemake@input[[2]]))
omega3 <- as.matrix(read.table(snakemake@input[[3]]))

# set up pdf output for the following plots
# all plots should be saved in the same pdf file
pdf(snakemake@output[[1]], width=8, height=8)
plot(omega1,omega2) ; abline(a=0,b=1) 
plot(omega1,omega3) ; abline(a=0,b=1)
plot(omega2,omega3) ; abline(a=0,b=1)
dev.off()

# compute the distance between the omega estimates
d12 = fmd.dist(omega1,omega2)
d13 = fmd.dist(omega1,omega3)
d23 = fmd.dist(omega2,omega3)

# save the distance in a file
write.table(data.frame(d12=d12, d13=d13, d23=d23),
            snakemake@output[[2]],
            sep=",", row.names=FALSE, quote=FALSE)

# STD IS model diagnostics
snp.res=read.table(snakemake@input[[1]],h=T) 

pdf(snakemake@output[[1]],width=10,height=10)
layout(matrix(1:3,3,1)) 
plot(snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)") 
plot(snp.res$eBPis,xlab="SNP",ylab="eBPis") 
plot(snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

# Contrast analysis diagnostics
# As authors describe in the manual:
# The population ecotype being a binary trait (either “crab” or “wave”), ## This refers to their example project!!! 
# one may rely on the C2 statistic (Olazcuaga et al., 2020) to identify SNPs 
# associated with the population ecotype rather than relying on the 
# (parametric) models used to estimate Bayes Factor.
ecotype.bf=snp.res$BF.dB. 
ecotype.C2=read.table(snakemake@input[[2]],h=T) 

#check the behavior of the p-values associated to the C2 
pdf(snakemake@output[[2]],width=10,height=10)
hist(10**(-1*ecotype.C2$log10.1.pval.),freq=F,breaks=50) 
abline(h=1) 
plot(ecotype.bf,ecotype.C2$log10.1.pval., 
    xlab="BF",ylab="C2 p-value (-log10 scale)") 
abline(h=3,lty=2) #0.001 p--value theshold 
abline(v=20,lty=2) #BF threshold for decisive evidence (according to Jeffreys’ rule)
dev.off()
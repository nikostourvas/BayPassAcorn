# STD IS model diagnostics
snp.res=read.table(snakemake@input[[1]],h=T)
envfactor_names = read.table(snakemake@input[[2]],h=F)

# Produce the diagnostic graphs suggested by Baypass with ggplot2, but additionally we will facet by COVARIABLE
library(ggplot2)

# each facet should get its name from the envfactor_names file
snp.res$COVARIABLE = envfactor_names[snp.res$COVARIABLE,]

# Manhattan plot of the Bayes factors
snp.res$COVARIABLE = factor(snp.res$COVARIABLE, levels = unique(snp.res$COVARIABLE))
plot_bf = ggplot(snp.res, aes(x = MRK, y = BF.dB.)) + 
  geom_point() +
  geom_point(data = snp.res[snp.res$BF.dB. > 20, ], aes(x = MRK, y = BF.dB.), color = "red") + 
  facet_wrap(~COVARIABLE, ncol = 1) +
  xlab("SNP") + ylab("BFis (in dB)")
ggsave(snakemake@output[[1]], plot = plot_bf, width = 10, height = 49, units = "in", dpi = 300)

# Manhattan plot of the eBPis
plot_ebpis = ggplot(snp.res, aes(x = MRK, y = eBPis)) +
    geom_point() +
    facet_wrap(~COVARIABLE, ncol = 1) +
    xlab("SNP") + ylab("eBPis")
ggsave(snakemake@output[[2]], plot = plot_ebpis, width = 10, height = 49, units = "in", dpi = 300)

# Manhattan plot of the Beta_is
plot_beta_is = ggplot(snp.res, aes(x = MRK, y = Beta_is)) +
    geom_point() +
    facet_wrap(~COVARIABLE, ncol = 1) +
    xlab("SNP") + ylab(expression(beta ~ "coefficient"))
ggsave(snakemake@output[[3]], plot = plot_beta_is, width = 10, height = 49, units = "in", dpi = 300)


# Genetic offset and RONA
# go.is = compute_genetic_offset(beta.coef=NULL,
#                                 regfile=snp.res,
#                                 covfile=snakemake@params[[1]], 
#                                 newenv='ADDRELEVANTFILE',refenv=NULL,
#                                 scalecov=TRUE,candidate.snp=NULL, 
#                                 compute.rona=TRUE)

# write.table(go.is$go,snakemake@output[[3]],quote=F,sep=",",row.names=F)
# write.table(go.is$BtBeigenvalues,snakemake@output[[4]],quote=F,sep=",",row.names=F)
# write.table(go.is$BtBeigenvectors,snakemake@output[[5]],quote=F,sep=",",row.names=F)
# write.table(go.is$covimp,snakemake@output[[6]],quote=F,sep=",",row.names=F)
# write.table(go.is$RONA,snakemake@output[[7]],quote=F,sep=",",row.names=F)
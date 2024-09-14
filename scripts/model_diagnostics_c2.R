# Contrast analysis diagnostics
# As authors describe in the manual:
# The population ecotype being a binary trait (either “crab” or “wave”), ## This refers to their example project!!! 
# one may rely on the C2 statistic (Olazcuaga et al., 2020) to identify SNPs 
# associated with the population ecotype rather than relying on the 
# (parametric) models used to estimate Bayes Factor.
snp.res=read.table(snakemake@input[[1]],h=T) 
ecotype.bf=snp.res$BF.dB. 
ecotype.C2=read.table(snakemake@input[[2]],h=T) 

# Generate the raster graphics
png("temp_hist.png", width = 10, height = 10, units = "in", res = 300)
hist(10**(-1*ecotype.C2$log10.1.pval.),freq=F,breaks=50) 
abline(h=1) 
dev.off()

png("temp_scatter.png", width = 10, height = 10, units = "in", res = 300)
plot(ecotype.bf,ecotype.C2$log10.1.pval., 
    xlab="BF",ylab="C2 p-value (-log10 scale)") 
abline(h=3,lty=2) #0.001 p--value theshold 
abline(v=20,lty=2) #BF threshold for decisive evidence (according to Jeffreys’ rule)
dev.off()

# Read the raster graphics back into R
library(png)
hist_img = readPNG("temp_hist.png")
scatter_img = readPNG("temp_scatter.png")

# Embed the raster graphics into the PDF
pdf(snakemake@output[[1]], width = 10, height = 10)
plot.new()
rasterImage(hist_img, 0, 0, 1, 1)
plot.new()
rasterImage(scatter_img, 0, 0, 1, 1)
dev.off()

# Clean up temporary files
file.remove("temp_hist.png", "temp_scatter.png")


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
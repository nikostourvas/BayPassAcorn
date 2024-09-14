# STD IS model diagnostics
snp.res=read.table(snakemake@input[[1]],h=T) 

png(snakemake@output[[1]],width=10,height=10, units="in", res=300)
layout(matrix(1:3,3,1)) 
plot(snp.res$BF.dB.,xlab="SNP",ylab="BFis (in dB)") 
plot(snp.res$eBPis,xlab="SNP",ylab="eBPis") 
plot(snp.res$Beta_is,xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

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
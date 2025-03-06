# re-direct console output to a log file for use with snakemake on HPC
# Open log file connection
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt") # 'wt' opens the file for writing text
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

# Description: Concatenate results from BayPass output files and snp_det files
# Modified from concatenate_res.R to include the option of retrieving Spearman correlation coefficients

concatenate_res2 <- function(dir="./",anaprefix="ana",extension="",nsubsets=2,
                          snpdet_prefix="./detsnp.sub",retrieve_pi_xtx=TRUE,retrieve_bfis=TRUE,retrieve_spearman=FALSE){
  #extension should be the same for snpdet files and baypass output files (i.e., all files should be compressed the same way or not)
  #snp_det_prefix should include the path (since it may be different than the directory containing BayPass output files)
  require(data.table)
  corepref=""
  if(nchar(anaprefix)>0){corepref=paste0(anaprefix,"_")}
  if(nsubsets<1){stop("Check nsubsets: at least one subset is needed\n")}
  
  for(i in 1:nsubsets){
    cat("Processing run",i,"out of",nsubsets,"\n")
    tmp.snpdet=fread((paste0(snpdet_prefix,i,extension)),data.table = F)[,1:2]
    colnames(tmp.snpdet)=c("CHR","POS")
    tmp.nsnps=nrow(tmp.snpdet)
    if(retrieve_pi_xtx){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_pi_xtx.out",extension),data.table=F)[,c("M_P","M_XtX","XtXst")]
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_bfis){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_betai_reg.out",extension),data.table=F)$"BF(dB)"
      tmp=matrix(as.numeric(tmp),tmp.nsnps)      
      colnames(tmp)=paste0("BFis_cov_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(retrieve_spearman){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_betai_reg.out",extension),data.table=F)$"M_Spearman"
      tmp=matrix(as.numeric(tmp),tmp.nsnps)      
      colnames(tmp)=paste0("M_Spearman",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(i==1){res=tmp.snpdet}else{res=rbind(res,tmp.snpdet)}
  }
  res=res[order(res$CHR,res$POS),]
  return(res)
}

library(ggplot2)
library(data.table)

all.res = concatenate_res2(
                        anaprefix = snakemake@params[[1]],
                        nsubsets=snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]], 
                        retrieve_pi_xtx=TRUE,
                        retrieve_bfis=TRUE,
                        retrieve_spearman=snakemake@params[[4]])

fwrite(all.res, snakemake@output[[1]], sep = ",", quote = FALSE)

print("Baypass results concatenated")

# Generate Manhattan plots per covariable
    # Produce the diagnostic graphs suggested by Baypass with ggplot
    all.res$CHR = as.factor(gsub("Qrob_Chr", "", all.res$CHR)) # remove prefix from chr names
    all.res$CHR = as.factor(gsub("H2.3_", "H2.3.", all.res$CHR)) # Fix for contig names
    all.res$CHR = as.factor(gsub("^0", "", all.res$CHR)) # replace 01, 02 etc with 1, 2 etc
    all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))

    # Manhattan plot of the Bayes factor for each covariable
    envfactor_names = read.table(snakemake@input[[1]],h=F)

    # Keep only the BF columns and throw away the Spearman columns
    all.res = all.res[ ,1:(5 + nrow(envfactor_names))] 

    # Produce the diagnostic graphs suggested by Baypass with ggplot2, but additionally we will facet by COVARIABLE
    # Ensure the number of columns in all.res matches the number of names in envfactor_names$V1
    if (ncol(all.res) != (5 + nrow(envfactor_names))) {
        stop("Mismatch between number of columns in all.res and number of names in envfactor_names")
    }

    # each facet should get its name from the envfactor_names file
    colnames(all.res) = c("CHR","POS","M_P","M_XtX", "XtXst", envfactor_names$V1)

    # generate tidy data set for plotting
    all.res = reshape2::melt(all.res, id.vars=c("CHR","POS","M_P","M_XtX", "XtXst"), 
                variable.name="COVARIABLE", value.name="BF.dB.")

    # Optional step
    # Remove SNPs with a Bayes factor lower than -10 to make plotting faster
    all.res = all.res[all.res$BF.dB. > 5,]

    # Manhattan plot of the Bayes factors
    p = ggplot(all.res,aes(x=POS, y=BF.dB., color=CHR)) + 
    geom_point(alpha=1, size=0.8) +
    geom_point(data=all.res[all.res$BF.dB.>20,], aes(x=POS, y=BF.dB.), color="red", size=.8) +
    facet_grid(COVARIABLE~CHR, space = "free_x", scales = "free_x") +
    scale_color_manual(values = rep(c("black", "grey"), 2000)) +
    theme_bw() +
    theme(legend.position = "none",
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.spacing.x = unit(0, "lines"),
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.minor.y = element_blank()) +
    xlab("Chromosome")+ ylab("BFis (in dB)")

    ggsave(snakemake@output[[2]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

print("Manhattan plots completed")

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection
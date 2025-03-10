library(ggplot2)
library(data.table)

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

all.res = concatenate_res2(
                        anaprefix = snakemake@params[[1]],
                        nsubsets=snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]], 
                        retrieve_pi_xtx=TRUE,
                        retrieve_bfis=TRUE,
                        retrieve_spearman=snakemake@params[[4]])


# Data preprocessing (edit this part to match the data structure of your experiment)
all.res$CHR = gsub("Qrob_Chr", "", all.res$CHR) # remove prefix from chr names
all.res$CHR = gsub("H2.3_", "H2.3.", all.res$CHR) # Fix for contig names
all.res$CHR = gsub("^0", "", all.res$CHR) # replace 01, 02 etc with 1, 2 etc
all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))

# Read environmental factor names
envfactor_names = read.table(snakemake@input[[1]], h=F)
print(paste("Found", nrow(envfactor_names), "environmental factors"))

# Check if column count matches before renaming
expected_cols = 5 + 2 * nrow(envfactor_names)  # CHR, POS, M_P, M_XtX, XtXst + BF and spearman columns
if(ncol(all.res) != expected_cols) {
  warning(paste("Expected", expected_cols, "columns but found", ncol(all.res)))
  # If column counts don't match, print column names for debugging
  print("Current column names:")
  print(colnames(all.res))
}

# Column renaming with explicit length check
new_colnames = c("CHR", "POS", "M_P", "M_XtX", "XtXst")
bf_names = paste0(envfactor_names$V1, "_BF")
spearman_names = paste0(envfactor_names$V1, "_spearman")
new_colnames = c(new_colnames, bf_names, spearman_names)

if(length(new_colnames) == ncol(all.res)) {
  colnames(all.res) = new_colnames
} else {
  warning("Column count mismatch. Using existing column names.")
}                        

# export BF
fwrite(all.res[, c("CHR", "POS", "M_P", "M_XtX", "XtXst", bf_names)],
       snakemake@output[[1]], sep = ",", quote = FALSE)

# generate absolute values of spearman rho and export
spearman_cols = all.res[, spearman_names]
spearman_cols = abs(spearman_cols)
spearman_cols = cbind(all.res[, c("CHR", "POS")], spearman_cols)
fwrite(spearman_cols, snakemake@output[[2]], sep = ",", quote = FALSE)
rm(spearman_cols)

print("Baypass results concatenated")

# Generate Manhattan plots per covariable
     # Keep only the BF columns and throw away the Spearman columns
     all.res = all.res[ ,1:(5 + nrow(envfactor_names))] 

     # generate tidy data set for plotting
     all.res = reshape2::melt(all.res, id.vars=c("CHR","POS","M_P","M_XtX", "XtXst"), 
                 variable.name="COVARIABLE", value.name="BF")

     # Optional step
     # Remove SNPs with a Bayes factor lower than -10 to make plotting faster
     all.res = all.res[all.res$BF > 5,]

    # Manhattan plot of the Bayes factors
    p = ggplot(all.res,aes(x=POS, y=BF, color=CHR)) + 
    geom_point(alpha=1, size=0.8) +
    geom_point(data=all.res[all.res$BF>20,], aes(x=POS, y=BF), color="red", size=.8) +
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
    xlab("Chromosome")+ ylab("BFis (dB)")

    ggsave(snakemake@output[[3]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

print("Manhattan plots completed")

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection
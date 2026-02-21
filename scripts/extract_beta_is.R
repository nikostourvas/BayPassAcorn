library(data.table)

# re-direct console output to a log file for use with snakemake on HPC
# Open log file connection
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt") # 'wt' opens the file for writing text
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

# Description: Extract chromosome, position, and regression coefficient Beta_is
# from BayPass output files. Reuses the same file-reading logic as concatenate_res2.R.

scipen=999  # avoid scientific notation
options(scipen=scipen)

concatenate_res2 <- function(dir="./", anaprefix="ana", extension="", nsubsets=2,
                             snpdet_prefix="./detsnp.sub", retrieve_pi_xtx=TRUE,
                             retrieve_bfis=TRUE, retrieve_spearman=FALSE,
                             retrieve_betais=FALSE){
  # extension should be the same for snpdet files and baypass output files
  # snp_det_prefix should include the path (since it may be different than the directory containing BayPass output files)
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
    if(retrieve_betais){
      tmp=fread(paste0(dir,"/",corepref,i,"_summary_betai_reg.out",extension),data.table=F)$"Beta_is"
      tmp=matrix(as.numeric(tmp),tmp.nsnps)
      colnames(tmp)=paste0("Beta_is_cov_",1:ncol(tmp))
      tmp.snpdet=cbind(tmp.snpdet,tmp)
    }
    if(i==1){res=tmp.snpdet}else{res=rbind(res,tmp.snpdet)}
  }
  res=res[order(res$CHR,res$POS),]
  return(res)
}

all.res = concatenate_res2(
                        anaprefix = snakemake@params[[1]],
                        nsubsets = snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]],
                        retrieve_pi_xtx = FALSE,
                        retrieve_bfis = FALSE,
                        retrieve_spearman = FALSE,
                        retrieve_betais = TRUE)

# Data preprocessing (edit this part to match the data structure of your experiment)
all.res$CHR = gsub("Qrob_Chr", "", all.res$CHR) # remove prefix from chr names
all.res$CHR = gsub("H2.3_", "H2.3.", all.res$CHR) # Fix for contig names
all.res$CHR = gsub("^0", "", all.res$CHR) # replace 01, 02 etc with 1, 2 etc
all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))

# Read environmental factor names
envfactor_names = read.table(snakemake@input[[1]], h=F)
print(paste("Found", nrow(envfactor_names), "environmental factors"))

# Check if column count matches before renaming
expected_cols = 2 + nrow(envfactor_names)  # CHR, POS + one Beta_is column per env factor
if(ncol(all.res) != expected_cols) {
  warning(paste("Expected", expected_cols, "columns but found", ncol(all.res)))
  print("Current column names:")
  print(colnames(all.res))
}

# Column renaming
betais_names = paste0(envfactor_names$V1, "_Beta_is")
new_colnames = c("CHR", "POS", betais_names)

if(length(new_colnames) == ncol(all.res)) {
  colnames(all.res) = new_colnames
} else {
  warning("Column count mismatch. Using existing column names.")
}

# Export Beta_is
fwrite(all.res, snakemake@output[[1]], sep = ",", quote = FALSE)

print("Beta_is extraction completed")

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection

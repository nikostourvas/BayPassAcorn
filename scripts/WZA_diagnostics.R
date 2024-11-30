library(ggplot2)
library(data.table)
library(reshape2)
library(qvalue)
library(tidyr)

FDR_LEVEL=snakemake@params[[2]]

############################
# Manhattan plots
############################

envfactors = readLines(snakemake@input[[1]])
envfactors = factor(envfactors, levels=envfactors)
print(envfactors)

WZA_res = list()
for (i in envfactors){
  WZA_res[[i]] = fread(paste0(snakemake@params[[1]], i, "_WZA_output.csv"))
}

# FDR correction
# for some reason, while testing with a toy dataset, WZA produced few NAs
# Until we investigate further, I replace those NAs with 1
WZA_res = lapply(WZA_res, function(x) {x$Z_pVal[is.na(x$Z_pVal)] = 1; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$qvalue = qvalue(x$Z_pVal, fdr.level=FDR_LEVEL)$significant; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$qvalue0.001 = qvalue(x$Z_pVal, fdr.level=0.001)$significant; return(x)})
# export each table to a file
for (i in envfactors){
  write.table(WZA_res[[i]], file=paste0(snakemake@params[[1]], i, "_WZA_output_fdr.csv"), 
              sep=",", row.names=FALSE, quote=FALSE)
}

WZA_res = lapply(WZA_res, function(x) x[ ,c(1,2,6,7,8,9)]) # keep only the columns we need for plotting
WZA_res = reshape2::melt(WZA_res, id.vars = c("index", "SNPs", "POS"))
# Make the dataset wide. Spread the "variable" column into multiple columns
WZA_res = spread(WZA_res, variable, value)
WZA_res$qvalue = as.logical(WZA_res$qvalue)
WZA_res$qvalue0.001 = as.logical(WZA_res$qvalue0.001)

WZA_res$index = as.factor(gsub("Qrob_Chr", "", WZA_res$index)) # remove prefix from chr names
WZA_res$index = as.factor(gsub("Sc0000", "", WZA_res$index)) # remove prefix from chromosome names
# remove from rows of the column WZA_res$index the suffix "_window_" and anything after that
WZA_res$index = as.factor(gsub("_window_.*", "", WZA_res$index))

colnames(WZA_res) = c("CHR", "SNPs", "POS", "envfactor", "value", "qvalue", "qvalue0.001")

WZA_res$CHR = as.factor(gsub("^0", "", WZA_res$CHR)) # replace 01, 02 etc with 1, 2 etc
WZA_res$CHR = factor(WZA_res$CHR, levels=unique(WZA_res$CHR))

# add a column with the calculated -log10(p-value)
WZA_res$score = -log10(WZA_res$value)

p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$qvalue==TRUE,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$qvalue0.001==TRUE,], aes(x=POS, y=score), color="red", size=.8) +
  facet_grid(factor(envfactor, levels=envfactors)~CHR, space = "free_x", scales = "free") +
  scale_color_manual(values = rep(c("black", "grey"), 2000)) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank()) +
  xlab("chromosome")+ ylab("-log10(q-value)")

ggsave(snakemake@output[[1]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

############################
# P-value distribution plots
############################

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_", env, ".png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$value, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$value))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$value)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}
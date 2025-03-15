library(ggplot2)
library(data.table)
library(reshape2)
library(qvalue)
library(tidyr)
library(viridis)

# re-direct console output to a log file for use with snakemake on HPC
# Open log file connection
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt") # 'wt' opens the file for writing text
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

FDR_LEVEL=snakemake@params[[2]]

# Read the environmental factor names
envfactors = readLines(snakemake@input[[1]])
envfactors = factor(envfactors, levels=envfactors)
print(envfactors)

WZA_res_BF = list()
for (i in envfactors){
  WZA_res_BF[[i]] = fread(paste0(snakemake@params[[1]], i, "_BF_WZA_output.csv"))[,c("index","POS","Z","Z_pVal")]
  colnames(WZA_res_BF[[i]]) = c("index","POS","Z_BF","Z_BF_pVal")
}

WZA_res_spearman = list()
for (i in envfactors){
  WZA_res_spearman[[i]] = fread(paste0(snakemake@params[[1]], i, "_spearman_WZA_output.csv"))[,c("index","POS","Z","Z_pVal")]
  colnames(WZA_res_spearman[[i]]) = c("index","POS","Z_spearman","Z_spearman_pVal")
}


# Merge the BF and Spearman results
WZA_res = list()
for (i in envfactors){
  WZA_res[[i]] = merge(WZA_res_BF[[i]], WZA_res_spearman[[i]], by=c("index", "POS"), all=TRUE)
}

# Print the first few rows of the merged data for debugging
print(head(WZA_res[[1]]))


# Calculate the genomic inflation factor (GIF) based on all the Z-scores.
gif_BF = list()
for (i in envfactors){
  gif_BF[[i]] = median((WZA_res[[i]]$Z_BF)^2)/(qchisq(0.5, df = 1, lower.tail = FALSE)) 
}
# Save the GIFs to a file
gif_table = reshape2::melt(gif_BF)
write.table(gif_table, file=paste0(snakemake@params[[1]], "GIF_BF.csv"), sep=",", row.names=FALSE, quote=FALSE)

gif_spearman = list()
for (i in envfactors){
  gif_spearman[[i]] = median((WZA_res[[i]]$Z_spearman)^2)/(qchisq(0.5, df = 1, lower.tail = FALSE)) 
}
# Save the GIFs to a file
gif_table = reshape2::melt(gif_spearman)
write.table(gif_table, file=paste0(snakemake@params[[1]], "GIF_spearman.csv"), sep=",", row.names=FALSE, quote=FALSE)

# Calculate the gif corrected p-values and q-values
print("Calculating p-values")
# for some reason, while testing with a toy dataset, WZA produced few NAs
# Until we investigate further, I replace those NAs with 1
WZA_res = lapply(WZA_res, function(x) {x$Z_BF_pVal[is.na(x$Z_BF_pVal)] = 1; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$Z_spearman_pVal[is.na(x$Z_spearman_pVal)] = 1; return(x)})
for (i in envfactors){
  WZA_res[[i]]$Z_BF_pVal_gif_adj = pchisq(WZA_res[[i]]$Z_BF^2/gif_BF[[i]], df = 1, lower.tail = FALSE)
  WZA_res[[i]]$Z_spearman_pVal_gif_adj = pchisq(WZA_res[[i]]$Z_spearman^2/gif_spearman[[i]], df = 1, lower.tail = FALSE)
}

print("Calculating q-values")
# FDR correction
WZA_res = lapply(WZA_res, function(x) {x$BF_qvalue = qvalue(x$Z_BF_pVal)$qvalues; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$spearman_qvalue = qvalue(x$Z_spearman_pVal)$qvalues; return(x)})

WZA_res = lapply(WZA_res, function(x) {x$BF_qvalue_gif_adj = qvalue(x$Z_BF_pVal_gif_adj)$qvalues; return(x)})
WZA_res = lapply(WZA_res, function(x) {x$spearman_qvalue_gif_adj = qvalue(x$Z_spearman_pVal_gif_adj)$qvalues; return(x)})

# export each table to a file
print("print tables to files")
for (i in envfactors){
  fwrite(WZA_res[[i]], paste0(snakemake@params[[1]], i, "_WZA_output_fdr.csv"))
}


############################
# Find common significant windows between BF and Spearman
############################

WZA_res = lapply(WZA_res, function(x) x[ ,c(1,2,4,6,7,8,9,10,11,12)]) # keep only the columns we need for plotting
# Combine the results into a single data frame
WZA_res = reshape2::melt(WZA_res, id.vars = c("index", "POS"))
# Make the dataset wide. Spread the "variable" column into multiple columns
WZA_res = spread(WZA_res, variable, value)
print(head(WZA_res))

print("string manipulation")
# remove from rows of the column WZA_res$index the suffix "_window_" and anything after that
WZA_res$index = as.factor(gsub("_window_.*", "", WZA_res$index))

colnames(WZA_res) = c("CHR", "POS", "envfactor", 
                      "BF_pvalue", "spearman_pvalue", "BF_gif_cor_pvalue", "spearman_gif_cor_pvalue", 
                      "BF_qvalue", "spearman_qvalue", "BF_qvalue_gif_adj", "spearman_qvalue_gif_adj")

WZA_res$CHR = factor(WZA_res$CHR, levels=unique(WZA_res$CHR))

print("Finding common significant windows between BF and Spearman")
# keep only rows where both BF and Spearman are significant (BF_qvalue_gif_adj < 0.05 and spearman_qvalue_gif_adj < 0.05)
WZA_sig0.1 = WZA_res[WZA_res$BF_qvalue_gif_adj < 0.1 & WZA_res$spearman_qvalue_gif_adj < 0.1,]
fwrite(WZA_sig0.1, snakemake@output[[5]])

WZA_sig0.05 = WZA_res[WZA_res$BF_qvalue_gif_adj < 0.05 & WZA_res$spearman_qvalue_gif_adj < 0.05,]
fwrite(WZA_sig0.05, snakemake@output[[6]])

WZA_sig0.01 = WZA_res[WZA_res$BF_qvalue_gif_adj < 0.01 & WZA_res$spearman_qvalue_gif_adj < 0.01,]
fwrite(WZA_sig0.01, snakemake@output[[7]])

WZA_sig0.001 = WZA_res[WZA_res$BF_qvalue_gif_adj < 0.001 & WZA_res$spearman_qvalue_gif_adj < 0.001,]
fwrite(WZA_sig0.001, snakemake@output[[8]])

WZA_p0.1 = WZA_res[WZA_res$BF_pvalue < 0.1 & WZA_res$spearman_pvalue < 0.1,]
fwrite(WZA_p0.1, snakemake@output[[9]])

WZA_p0.05 = WZA_res[WZA_res$BF_pvalue < 0.05 & WZA_res$spearman_pvalue < 0.05,]
fwrite(WZA_p0.05, snakemake@output[[10]])

WZA_p0.01 = WZA_res[WZA_res$BF_pvalue < 0.01 & WZA_res$spearman_pvalue < 0.01,]
fwrite(WZA_p0.01, snakemake@output[[11]])

WZA_p0.001 = WZA_res[WZA_res$BF_pvalue < 0.001 & WZA_res$spearman_pvalue < 0.001,]
fwrite(WZA_p0.001, snakemake@output[[12]])

############################
# Manhattan plots
############################
print("Creating Manhattan plots")

# Sort the CHR column properly
# First remove prefixes from contig names
WZA_res$CHR = gsub("^Qrob_H2.3.Sc00000", "", WZA_res$CHR)
WZA_res$CHR = gsub("^Qrob_H2.3.Sc0000", "", WZA_res$CHR)
# Then convert to numeric
WZA_res$CHR = as.integer(WZA_res$CHR)
# Then convert back to factor with properly sorted levels
# If I don't do this, scale_color_manual will not work properly
WZA_res$CHR = factor(WZA_res$CHR, levels=as.character(sort(unique(WZA_res$CHR))))
print(unique(WZA_res$CHR))

# add a column with the calculated -log10(p-value)
WZA_res$score = -log10(WZA_res$BF_gif_cor_pvalue)
p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$BF_qvalue_gif_adj < 0.05,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$BF_qvalue_gif_adj < 0.001,], aes(x=POS, y=score), color="red", size=.8) +
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
  xlab("chromosome")+ ylab("-log10(p-value)")
ggsave(snakemake@output[[1]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

WZA_res$score = -log10(WZA_res$BF_pvalue)
p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$BF_qvalue < 0.05,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$BF_qvalue < 0.001,], aes(x=POS, y=score), color="red", size=.8) +
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
  xlab("chromosome")+ ylab("-log10(p-value)")
ggsave(snakemake@output[[2]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

WZA_res$score = -log10(WZA_res$spearman_gif_cor_pvalue)
p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$spearman_qvalue_gif_adj < 0.05,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$spearman_qvalue_gif_adj < 0.001,], aes(x=POS, y=score), color="red", size=.8) +
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
  xlab("chromosome")+ ylab("-log10(p-value)")
ggsave(snakemake@output[[3]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)

WZA_res$score = -log10(WZA_res$spearman_pvalue)
p <- ggplot(WZA_res, aes(x=POS, y=score, color=CHR)) + 
  geom_point(alpha=1, size=0.8) +
  geom_point(data=WZA_res[WZA_res$spearman_qvalue < 0.05,], aes(x=POS, y=score), color="darkred", size=.8) +
  geom_point(data=WZA_res[WZA_res$spearman_qvalue < 0.001,], aes(x=POS, y=score), color="red", size=.8) +
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
  xlab("chromosome")+ ylab("-log10(p-value)")
ggsave(snakemake@output[[4]],width=10,height=60, units="in", dpi=300, limitsize=FALSE)


############################
# P-value distribution plots
############################
print("Creating P-value distribution plots")

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_", env, "_BF.png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$BF_pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$BF_pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$BF_pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_(GIF_corrected)_", env, "_BF.png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$BF_gif_cor_pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$BF_gif_cor_pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$BF_gif_cor_pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_(Q-value)_", env, "_spearman.png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))

  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$spearman_pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$spearman_pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$spearman_pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}

# Loop over each envfactor
for (env in envfactors) {
  # Set up the file name for the PNG
  png_filename <- paste0(snakemake@params[[1]], "Pvalue_and_QQ_Plots_(GIF_corrected)_", env, "_spearman.png")
  
  # Open the PNG device
  png(filename = png_filename, width = 1050, height = 750, units = "px")
  
  # Set up the plotting area
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 1))
  
  # Subset data for the current envfactor
  data_subset <- subset(WZA_res, envfactor == env)
  
  # Plot the p-value distribution
  hist(data_subset$spearman_gif_cor_pvalue, col = "grey", main = paste("P-value Distribution -", env), xlab = "P-values")
  
  # Plot the QQ-plot
  observed <- -log10(sort(data_subset$spearman_gif_cor_pvalue))
  expected <- -log10(ppoints(length(observed)))
  plot(expected, observed, xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", main = paste("QQ-Plot -", env), pch = 19, cex = 1)
  abline(0, 1)
  
  # Add median p-value to the QQ-plot
  median_pvalue <- median(data_subset$spearman_gif_cor_pvalue)
  points(-log10(median_pvalue), -log10(median_pvalue), col = "red", pch = 17, cex = 2)
  legend("bottomright", legend = paste("Median P = ", round(median_pvalue, 4)), col = "red", pch = 17, cex = 1)
  
  # Close the PNG device
  dev.off()
}

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection
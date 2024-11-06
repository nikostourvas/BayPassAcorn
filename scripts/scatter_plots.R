# Produce scatter plots of genotype frequencies vs environmental values
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(forcats)

# Load the concatenated Baypass results
all.res = fread(snakemake@input[[1]])
all.res$CHR = as.factor(gsub("Qrob_Chr", "", all.res$CHR)) # remove prefix from chr names
all.res$CHR = as.factor(gsub("^0", "", all.res$CHR)) # replace 01, 02 etc with 1, 2 etc
all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))
all.res$CHR = as.numeric(all.res$CHR)

envfactor_names = read.table(snakemake@input[[3]],h=F)
# each facet should get its name from the envfactor_names file
colnames(all.res) = c("CHR","POS","M_P","M_XtX", "XtXst", envfactor_names$V1)

# generate tidy data set for plotting
all.res = melt(all.res, id.vars=c("CHR","POS","M_P","M_XtX", "XtXst"), 
               variable.name="COVARIABLE", value.name="BF.dB.")

# Extract significant SNPs
all.res$CHR = as.numeric(all.res$CHR)
sig.snps = all.res[all.res$BF.dB. > snakemake@params[[1]],]

# Extract allele frequencies
afs = fread(snakemake@input[[2]], header=TRUE)
# Remove the "Qrob_Chr" prefix from the chrom_pos column
afs[, chrom_pos := gsub("Qrob_Chr", "", chrom_pos)]
# Split the chrom_pos column by the underscore ('_')
afs[, c("CHR", "POS") := tstrsplit(chrom_pos, "_", keep = c(1, 2))]
# Convert the CHR and POS columns to numeric
afs[, CHR := as.numeric(CHR)]
afs[, POS := as.numeric(POS)]
# make it tidy
afs_tidy = reshape2::melt(afs, id.vars=c("chrom_pos", "CHR","POS"), variable.name="POP", value.name="AF")

# Keep only the rows of afs that correspond to the significant SNPs based on the CHR and POS columns
afs_tidy = merge(afs_tidy, sig.snps, by=c("CHR","POS"))

# Load environmental data
env_factor_names = read.table(snakemake@input[[3]],h=F)
env_factor_data = read.table(snakemake@input[[4]],h=F)
env_factors = cbind(env_factor_names, env_factor_data)
rm(env_factor_names, env_factor_data)
colnames(env_factors) = c("COVARIABLE", colnames(afs)[2:(ncol(afs)-2)])
# make it tidy
env_factors = reshape2::melt(env_factors, id.vars="COVARIABLE", variable.name="POP", value.name="env_value")

# Merge environmental data with Allele Frequencies of significant SNPs
fulldata = merge(afs_tidy, env_factors, by=c("COVARIABLE", "POP"))

# Plot scatter plots of allele frequencies vs environmental values for each covariable and each chrom_pos
for (covar in unique(fulldata$COVARIABLE)) { 
    covariable_dat = fulldata[fulldata$COVARIABLE == covar,]
    
    # Check if covariable_dat is empty
    if (nrow(covariable_dat) == 0) {
        next
    }
    
    # Reorder `facet_var` based on `order_var` (e.g., using mean as criterion)
    covariable_dat$chrom_pos <- fct_reorder(covariable_dat$chrom_pos, 
                                            covariable_dat$BF.dB., 
                                            .fun = mean, .desc = TRUE)

    # Systematically subsample the SNP list to retain only 150 plots or less
    max_plots = 150
    if (length(unique(covariable_dat$chrom_pos)) > max_plots) {
        step = ceiling(length(unique(covariable_dat$chrom_pos)) / max_plots)
        chrom_pos_to_keep = unique(covariable_dat$chrom_pos)[seq(1, length(unique(covariable_dat$chrom_pos)), step)]
        covariable_dat = covariable_dat[covariable_dat$chrom_pos %in% chrom_pos_to_keep,]
    }
    
    # Summarize data to create a unique label for each facet
    label_dat <- covariable_dat %>%
        group_by(chrom_pos) %>%
        reframe(label_text = paste("M_P:", round(M_P, 2),"\n",
                                    "M_XtX:", round(M_XtX, 2),"\n",
                                    "XtXst:", round(XtXst, 2),"\n",
                                    "BF.dB:", round(BF.dB., 2)))
    
    p = ggplot(covariable_dat, aes(x=env_value, y=AF)) +
        geom_point(alpha=0.5, color="blue") +
        scale_y_continuous(limits = c(0, 1)) +
        facet_wrap(~chrom_pos, ncol=2) +
        geom_text(data=label_dat, aes(label=label_text), 
                        x=max(covariable_dat$env_value), y=0.95, size=2,
                        hjust=1, vjust=1) +
        coord_cartesian(clip = "off") +
        theme_bw() +
        theme(legend.position = "none", plot.title = element_text(size = 8)) +
        xlab(unique(covariable_dat$COVARIABLE)) + ylab("Genotype frequency")
    
    # Dynamic Scaling for Large Numbers of Plots 
    # Determine the number of unique facets
    num_facets <- length(unique(covariable_dat$chrom_pos))
    
    # Calculate grid dimensions
    cols <- 2  # Fixed number of columns, for example
    rows <- ceiling(num_facets / cols)
    
    # Adjust dimensions based on desired width and height per panel
    panel_width <- 3   # width of each panel in inches
    panel_height <- 3  # height of each panel in inches
    
    total_width <- panel_width * cols
    total_height <- panel_height * rows
    
    # Save the plot with calculated dimensions
    plot_filename <- paste0(snakemake@params[[2]], covar, "_scatterplots.pdf")
    ggsave(plot_filename, p, width = total_width, height = total_height, units = "in", limitsize = FALSE)
}

print("Scatter plots completed")
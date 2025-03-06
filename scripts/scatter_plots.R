# Produce scatter plots of genotype frequencies vs environmental values
library(ggplot2)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(forcats)

# re-direct console output to a log file for use with snakemake on HPC
# Open log file connection
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt") # 'wt' opens the file for writing text
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

############################################
# Load the concatenated Baypass results
############################################

all.res = fread(snakemake@input[[1]])
print(paste("Loaded data with", ncol(all.res), "columns and", nrow(all.res), "rows"))

# Data preprocessing
all.res$CHR = as.factor(gsub("Qrob_Chr", "", all.res$CHR))
all.res$CHR = as.factor(gsub("H2.3_", "H2.3.", all.res$CHR)) # Fix for contig names
all.res$CHR = as.factor(gsub("^0", "", all.res$CHR))
all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))
all.res$CHR = as.numeric(all.res$CHR)

# Read environmental factor names
envfactor_names = read.table(snakemake@input[[3]], h=F)
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

# Check for BF columns
bf_cols = grep("_BF$", names(all.res), value=TRUE)
if(length(bf_cols) == 0) {
  stop("No BF columns found. Check column naming.")
}

# Generate tidy dataset
print(paste("Found", length(bf_cols), "BF columns"))
all.res_subset_BF = all.res[, c("CHR", "POS", "M_P", "M_XtX", "XtXst", bf_cols), with=FALSE]
all.res_subset_spearman = all.res[, c("CHR", "POS", spearman_names), with=FALSE]

all.res_subset_BF = melt(all.res_subset_BF, 
               id.vars=c("CHR", "POS", "M_P", "M_XtX", "XtXst"),
               variable.name="COVARIABLE", value.name="BF.dB.")
# remove the "_BF" suffix from the COVARIABLE column
all.res_subset_BF$COVARIABLE = gsub("_BF", "", all.res_subset_BF$COVARIABLE)

all.res_subset_spearman = melt(all.res_subset_spearman, 
               id.vars=c("CHR", "POS"),
               variable.name="COVARIABLE", value.name="spearman_rho")

all.res = cbind(all.res_subset_BF, all.res_subset_spearman$spearman_rho)
colnames(all.res) = c("CHR", "POS", "M_P", "M_XtX", "XtXst", "COVARIABLE", "BF.dB.", "spearman_rho")

# remove redundant datasets
rm(all.res_subset_BF, all.res_subset_spearman)

# remoce low signal SNPs to lighten processing load
all.res = all.res[all.res$BF.dB. > 5,]

############################################
# Extract significant SNPs
############################################
all.res$CHR = as.numeric(all.res$CHR)
sig.snps = all.res[all.res$BF.dB. > snakemake@params[[1]] 
                   & abs(all.res$spearman_rho) >= snakemake@params[[2]], ]
print(paste("Found", nrow(sig.snps), "significant SNPs"))

# Exit early if no significant SNPs
if(nrow(sig.snps) == 0) {
  message("No SNPs exceed the BF threshold. Creating empty output files.")
  # Create empty output files to satisfy Snakemake
  for(covar in unique(envfactor_names$V1)) {
    plot_filename <- paste0(snakemake@params[[3]], covar, "_scatterplots.pdf")
    pdf(plot_filename, width=8, height=6)
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    text(1, 1, "No significant SNPs found", cex=1.5)
    dev.off()
  }
}

# Plot BF vs Spearman rho
p = ggplot(all.res, aes(x=spearman_rho, y=BF.dB.)) +
  geom_point(alpha=0.5, color="darkblue") +
  geom_point(alpha=0.5, data=sig.snps, aes(x=spearman_rho, y=BF.dB.), color="darkred") +
  geom_hline(yintercept = snakemake@params[[1]], linetype="dashed") +
  geom_vline(xintercept = snakemake@params[[2]], linetype="dashed") +
  geom_vline(xintercept = -snakemake@params[[2]], linetype="dashed") +
  xlab("Spearman rho") + ylab("BF (db)") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(size = 8))

ggsave(snakemake@output[[1]], p, width = 8, height = 6, units = "in")

############################################
# Load allele frequencies
############################################
afs = fread(snakemake@input[[2]], header=TRUE)
# Remove the "Qrob_Chr" prefix from the chrom_pos column
afs[, chrom_pos := gsub("Qrob_Chr", "", chrom_pos)]
afs[, chrom_pos := gsub("H2.3_", "H2.3.", chrom_pos)] # Fix for contig names
# Split the chrom_pos column by the underscore ('_')
afs[, c("CHR", "POS") := tstrsplit(chrom_pos, "_", keep = c(1, 2))]

print("convert to numeric")
# Convert the CHR and POS columns to numeric
afs[, CHR := as.numeric(CHR)]
afs[, POS := as.numeric(POS)]

print("tidy data")
# make it tidy
afs_tidy = reshape2::melt(afs, id.vars=c("chrom_pos", "CHR","POS"), variable.name="POP", value.name="AF")

print("filter (merge) data")
# Keep only the rows of afs that correspond to the significant SNPs based on the CHR and POS columns
afs_tidy = merge(afs_tidy, sig.snps, by=c("CHR","POS"))

print("AF data loaded")

############################################
# Load environmental data
############################################
env_factor_names = read.table(snakemake@input[[3]],h=F)
env_factor_data = read.table(snakemake@input[[4]],h=F)
env_factors = cbind(env_factor_names, env_factor_data)
rm(env_factor_names, env_factor_data)
colnames(env_factors) = c("COVARIABLE", colnames(afs)[2:(ncol(afs)-2)])
# make it tidy
env_factors = reshape2::melt(env_factors, id.vars="COVARIABLE", variable.name="POP", value.name="env_value")

print("Environmental data loaded")

# Merge environmental data with Allele Frequencies of significant SNPs
fulldata = merge(afs_tidy, env_factors, by=c("COVARIABLE", "POP"))

print("Data merged")

# output the fulldata table
fwrite(fulldata, snakemake@output[[2]])

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
                                    "BF.dB:", round(BF.dB., 2),"\n",
                                    "rho:", round(spearman_rho, 2)))
    
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
    plot_filename <- paste0(snakemake@params[[3]], covar, "_scatterplots.pdf")
    ggsave(plot_filename, p, width = total_width, height = total_height, units = "in", limitsize = FALSE)
}

print("Scatter plots completed")

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection
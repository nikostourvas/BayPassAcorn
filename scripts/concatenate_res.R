source("scripts/baypass_utils.R")

all.res = concatenate_res(
                        anaprefix = snakemake@params[[1]],
                        nsubsets=snakemake@params[[2]],
                        snpdet_prefix = snakemake@params[[3]], 
                        retrieve_pi_xtx=TRUE,
                        retrieve_bfis=TRUE,
                        retrieve_c2=snakemake@params[[4]])

write.table(all.res, snakemake@output[[1]], sep = ",", 
            row.names = FALSE, quote = FALSE)

# if snakemake@params[[4]] is TRUE, script is finished. Otherwise, generate Manhattan plots per covariable
if (!snakemake@params[[4]]) {
    # Produce the diagnostic graphs suggested by Baypass with ggplot
    all.res$CHR = as.factor(gsub("Qrob_Chr", "", all.res$CHR)) # remove prefix from chr names
    all.res$CHR = as.factor(gsub("^0", "", all.res$CHR)) # replace 01, 02 etc with 1, 2 etc
    all.res$CHR = factor(all.res$CHR, levels=unique(all.res$CHR))

    # Manhattan plot of the Bayes factor for each covariable
    envfactor_names = read.table(snakemake@input[[1]],h=F)

    # Produce the diagnostic graphs suggested by Baypass with ggplot2, but additionally we will facet by COVARIABLE
    library(ggplot2)
    library(data.table)

    # Ensure the number of columns in all.res matches the number of names in envfactor_names$V1
    if (ncol(all.res) != (5 + nrow(envfactor_names))) {
        stop("Mismatch between number of columns in all.res and number of names in envfactor_names")
    }

    # each facet should get its name from the envfactor_names file
    colnames(all.res) = c("CHR","POS","M_P","M_XtX", "XtXst", envfactor_names$V1)

    # generate tidy data set for plotting
    all.res = melt(all.res, id.vars=c("CHR","POS","M_P","M_XtX", "XtXst"), 
                variable.name="COVARIABLE", value.name="BF.dB.")

    # Optional step
    # Remove SNPs with a Bayes factor lower than -10 to make plotting faster
    all.res = all.res[all.res$BF.dB. > -10,]

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
}
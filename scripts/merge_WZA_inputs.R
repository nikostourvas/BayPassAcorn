# re-direct console output to a log file for use with snakemake on HPC
# Open log file connection
log_file <- snakemake@log[[1]]
log_conn <- file(log_file, open = "wt") # 'wt' opens the file for writing text
sink(log_conn, append = TRUE)
sink(log_conn, append = TRUE, type = "message")

library(GenomicRanges)
library(data.table)
library(dplyr)
library(rtracklayer)

# Read the input files
tmp1 <- fread(snakemake@input[[1]])
tmpMAF <- fread(snakemake@input[[2]])

# Extract chromosome prefix from the Identifier column in tmpMAF
# Example: "Qrob_Chr01_2087" -> extract just "Qrob_Chr" (without the number)
# We'll get the full format by removing the position suffix
example_chr <- sub("_[^_]+$", "", tmpMAF$Identifier[1])  # "Qrob_Chr01"
chr_prefix <- sub("[0-9]+$", "", example_chr)  # "Qrob_Chr"

# Normalize chromosome names in tmp1
# Only replace the specific pattern ".Sc" with "_Sc" to match tmpMAF format
# "Qrob_H2.3.Sc0000515" -> "Qrob_H2.3_Sc0000515"
tmp1$CHR_normalized <- gsub("\\.Sc", "_Sc", tmp1$CHR)

# Format chromosome names in tmp1 to match tmpMAF format
# For regular chromosomes: "1" -> "Qrob_Chr01"
# For scaffolds: "Qrob_H2.3_Sc0000515" -> "Qrob_H2.3_Sc0000515" (already formatted)
tmp1$CHR_formatted <- ifelse(
    grepl("^[0-9]+$", tmp1$CHR_normalized),  # If it's just a number
    paste0(chr_prefix, sprintf("%02d", as.integer(tmp1$CHR_normalized))),  # Format as chromosome
    tmp1$CHR_normalized  # Otherwise keep the scaffold name as is
)

# Format chromosome names in tmpMAF by extracting from Identifier
# Example: "Qrob_Chr01_2087" -> "Qrob_Chr01"
# Example: "Qrob_H2.3_Sc0000515_505014" -> "Qrob_H2.3_Sc0000515"
tmpMAF$CHR_formatted <- sub("_[^_]+$", "", tmpMAF$Identifier)

# Extract position from tmpMAF Identifier
# Example: "Qrob_Chr01_2087" -> 2087
tmpMAF$POS <- as.integer(sub(".*_", "", tmpMAF$Identifier))

# Create GRanges objects
gr_tmp1 <- GRanges(
    seqnames = tmp1$CHR_formatted,
    ranges = IRanges(start = tmp1$POS, width = 1)
)

gr_maf <- GRanges(
    seqnames = tmpMAF$CHR_formatted,
    ranges = IRanges(start = tmpMAF$POS, width = 1)
)

# Check if GFF file is provided and exclude regions
gff_file <- snakemake@params[["gff"]]
if (!is.null(gff_file) && gff_file != "" && file.exists(gff_file)) {
    cat("Loading regions from GFF file:", gff_file, "\n")
    
    # Import GFF file using rtracklayer
    # This automatically handles coordinate conversions and GFF format parsing
    gr_gff <- import(gff_file, format = "gff")
    
    cat("Total regions loaded from GFF:", length(gr_gff), "\n")
    
    # Find SNPs that overlap with GFF regions
    gff_overlaps_tmp1 <- findOverlaps(gr_tmp1, gr_gff)
    gff_overlaps_maf <- findOverlaps(gr_maf, gr_gff)
    
    # Get indices of SNPs NOT overlapping with GFF regions
    tmp1_keep <- setdiff(seq_along(gr_tmp1), queryHits(gff_overlaps_tmp1))
    maf_keep <- setdiff(seq_along(gr_maf), queryHits(gff_overlaps_maf))
    
    cat("SNPs in tmp1 overlapping GFF regions:", length(queryHits(gff_overlaps_tmp1)), "\n")
    cat("SNPs in MAF overlapping GFF regions:", length(queryHits(gff_overlaps_maf)), "\n")
    
    # Subset GRanges to exclude GFF regions
    gr_tmp1 <- gr_tmp1[tmp1_keep]
    gr_maf <- gr_maf[maf_keep]
    
    # Also subset the data.tables
    tmp1 <- tmp1[tmp1_keep, ]
    tmpMAF <- tmpMAF[maf_keep, ]
    
    cat("SNPs retained after GFF filtering - tmp1:", length(gr_tmp1), "\n")
    cat("SNPs retained after GFF filtering - MAF:", length(gr_maf), "\n")
} else {
    cat("No GFF file provided or file does not exist. Skipping filtering.\n")
}

# Debug - print first few unique chromosomes to verify formatting
cat("\nFirst few unique chromosomes in tmp1:\n")
print(head(unique(tmp1$CHR)))
cat("\nFirst few normalized chromosomes in tmp1:\n")
print(head(unique(tmp1$CHR_normalized)))
cat("\nFirst few formatted chromosomes in tmp1:\n")
print(head(unique(tmp1$CHR_formatted)))
cat("\nLast few formatted chromosomes in tmp1 (scaffolds):\n")
print(tail(unique(tmp1$CHR_formatted)))
cat("\nFirst few formatted chromosomes in tmpMAF:\n")
print(head(unique(tmpMAF$CHR_formatted)))
cat("\nLast few formatted chromosomes in tmpMAF (scaffolds):\n")
print(tail(unique(tmpMAF$CHR_formatted)))

# Find overlaps
overlaps <- findOverlaps(gr_tmp1, gr_maf, type = "equal")

# Create matched dataframes
tmp1_matched <- tmp1[queryHits(overlaps), ]
maf_matched <- tmpMAF$MAF[subjectHits(overlaps)]

# Combine the data (remove temporary columns used for matching)
tmp1_matched$CHR_normalized <- NULL
tmp1_matched$CHR_formatted <- NULL
result <- cbind(tmp1_matched, MAF = maf_matched)

# Write output
fwrite(result, snakemake@output[[1]], sep = ",")

# Report matching statistics
cat(sprintf("\nTotal SNPs in tmp1: %d\n", nrow(tmp1)))
cat(sprintf("Total SNPs in MAF table: %d\n", nrow(tmpMAF)))
cat(sprintf("Matched SNPs: %d\n", nrow(result)))
cat(sprintf("Unmatched SNPs from tmp1: %d\n", nrow(tmp1) - nrow(result)))

# Check if there are unmatched SNPs
if (nrow(result) < nrow(tmp1)) {
    warning(sprintf("%d SNPs from tmp1 could not be matched to MAF table", 
                    nrow(tmp1) - nrow(result)))
}

# Close the log file connection at the end of the script
sink() # Close standard output redirection
sink(type = "message") # Close message redirection
close(log_conn) # Close the file connection
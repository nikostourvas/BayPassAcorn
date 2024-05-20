#!/bin/bash

# Nikos Tourvas & Lars Littmann 
# 2023.11.14
# Take a VCF produced by VarScan and output a table containing the allele frequencies for all samples per position.

# declare variables
VCF=${1}

#Take the name of the input vcf without file type suffix.
OUTPUT=${VCF/.vcf.gz/}

#Extract the number of samples found in the vcf.
N_SAMPLES=$(bcftools query -l ${VCF} | wc -l)

# We extract the following information from the VCF: 
# Chromosome, Position, Ref Allele, Alt Allele, Ref Allele read depth (per pool), Alt Allele read depth (per pool)
# bcftools query arguments (extract data from VCF)
# -H: include header
# -f: extract FORMAT columns

# After the table is produced, it is fed to awk for transformations. These transformations are:
# 1.Fix table header: (a) remove square brackets and replace first four column names 
#                     (b) replace the suffix ":RD" with ".ref.cnt" in columns which signify ref allele counts
#                     (c) replace the suffix ":AD" with ".alt.cnt" in columns which signify alt allele counts
# 2.Replace periods (.) which signify missing data with zeros (0) in the allele count columns

bcftools query -Hf '%CHROM\t%POS\t%REF\t%ALT[\t%RD][\t%AD]\n' ${VCF}\
   | awk 'BEGIN{OFS=FS="\t"} 
        NR==1 {
            $1="chrom"; $2="pos"; $3="ref"; $4="alt";
            for (i=5; i<=NF; i++){
                gsub(/\[[^]]*\]/, "", $i); 
                gsub(/:RD$/, ".ref.cnt", $i); 
                gsub(/:AD$/, ".alt.cnt", $i)
            };
        print; next} 
        NR>1 {for (i=2; i<=NF; i++) 
            gsub(/\./, NA, $i); print
        }' \
| sed -e 's/\t\t/\tNA\t/g' > ${OUTPUT}_counts.txt
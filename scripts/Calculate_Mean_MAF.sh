#!/usr/bin/awk -f

# Description: Calculates the mean allele frequency and Mean Minor Allele Frequency (MAF) for each row.

BEGIN {
    FS = "\t"
    OFS = ","
}

NR == 1 {
    # Print headers for the output
    print "Identifier", "Mean_AF", "MAF"
}

NR > 1 {
    sum = 0
    count = 0
    for (i = 2; i <= NF; i++) {
        sum += $i
        count++
    }
    mean = sum / count
    maf = (mean <= 0.5) ? mean : 1 - mean
    # Format output to 3 decimals with tab separators
    printf "%s,%.3f,%.3f\n", $1, mean, maf
}

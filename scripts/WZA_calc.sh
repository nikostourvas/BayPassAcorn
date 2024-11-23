#!/bin/bash

# Nikolaos Tourvas
# 2021-09-29
# This script is a draft pipeline for running the WZA analysis based on BayPass output.
# Make_Genomic_Windows.sh separates the genome into windows of a specified size.
# Calculate_Mean_MAF.sh calculates the mean allele frequency and MAF for each SNP.

STATISTIC_FILE=$1 # Set the summary statistic to use for WZA
PREFIX=$(basename ${STATISTIC_FILE} .csv)
WINDOW_SIZE=$2 # Set the window size for WZA
COVARIABLE_STATISTIC=$3

awk -f scripts/Make_Genomic_Windows.sh -v window_size=${WINDOW_SIZE} results/${PREFIX}.csv > tmp1

awk -f scripts/Calculate_Mean_MAF.sh data/${PREFIX}_AlleleFrequencyTable.txt > tmp2

cut -d',' -f3 tmp2 > tmpMAF

paste -d',' tmp1 tmpMAF > results/${PREFIX}_WZA_input.csv

rm tmp1 tmp2 tmpMAF

python3 scripts/general_WZA_script.py \
    --correlations results/${PREFIX}_WZA_input.csv \
    --summary_stat BFis_${COVARIABLE_STATISTIC} \
    --window window_id \
    --MAF MAF \
    --sep "," \
    --retain POS \
    --output results/${PREFIX}_WZA_output.csv
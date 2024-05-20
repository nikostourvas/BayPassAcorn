#!/bin/bash

# Define the input file
INPUT=$1

# Check if the file exists
if [ ! -f ${INPUT} ]; then
    echo "File not found!"
    exit 1
fi

# Replace the string
sed 's/P01.....ACORN.BOKU.Pl..//g' ${INPUT}
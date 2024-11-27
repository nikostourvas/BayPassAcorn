#!/usr/bin/awk -f

#Assisgn SNPs to genomic windows of a user-defined size. By specifying the 
#window size at runtime, you can easily adjust your analysis without modifying 
#the script itself. 

BEGIN {
    FS = ",";  # Set the input field separator to tab
    OFS = ","; # Set the output field separator to tab

    # Validate the 'window_size' variable provided by the user
    if (!window_size || window_size !~ /^[0-9]+$/ || window_size <= 0) {
        # If 'window_size' is missing, not a number, or not positive, print an error and exit
        print "Error: Please provide a valid window size (positive integer) using -v window_size=VALUE" > "/dev/stderr";
        exit 1;
    }
}

# Process the header line (NR == 1)
NR == 1 {
    # Loop through all fields in the header to find 'CHR' and 'POS' columns
    for (i = 1; i <= NF; i++) {
        if ($i == "CHR") {
            CHR_col = i;  # Store the index of the 'CHR' column
        }
        if ($i == "POS") {
            POS_col = i;  # Store the index of the 'POS' column
        }
    }

    # Check if both 'CHR' and 'POS' columns were found
    if (!CHR_col || !POS_col) {
        # If not, print an error message and exit
        print "Error: Columns 'CHR' and 'POS' not found." > "/dev/stderr";
        exit 1;
    }

    # Print the header line, adding the new columns 'window_number' and 'window_id'
    print $0, "window_number", "window_id";
    next;  # Move to the next line
}

# Process the data lines (for NR > 1)
{
    # Extract the 'CHR' and 'POS' values using their column indices
    CHR = $CHR_col;
    POS = $POS_col;

    # Check if 'POS' is a numeric value
    if (POS ~ /^[0-9]+$/) {
        # Calculate the window number
        # Subtract 1 from 'POS' because POSitions start at 1
        # Divide by 'window_size' to get the window index
        # Use 'int' to get the integer part (floor division)
        # Add 1 to start window numbering from 1
        window_number = int((POS - 1) / window_size) + 1;
        # Zero-pad the window_number to 5 digits
        window_number_padded = sprintf("%05d", window_number);
    } else {
        # If 'POS' is not numeric, set 'window_number' to 'NA'
        window_number_padded = "NA";
    }

    # Construct the 'window_id' by combining 'CHR' and 'window_number'
    window_id = CHR "_window_" window_number_padded;

    # Print the original data line with the new window information
    print $0, window_number_padded, window_id;
}


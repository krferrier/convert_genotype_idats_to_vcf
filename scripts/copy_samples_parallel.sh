#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <n_threads> <valid_samples.txt> <idat-dir> <output-dir>"
    exit 1
fi

# Assign command line arguments to variables
n_threads="$1"
sample_file="$2"
raw_data_dir="$3"
destination_directory="$4"

# Check if the destination directory exists, create it if not
if [ ! -d "$destination_directory" ]; then
    mkdir -p "$destination_directory"
fi

# Function to copy .idat files for a given sample
copy_idats() {
    sample="$1"
    dest_dir="$2"
    source_idats="$3"

    find "$source_idats" -type f -name "${sample}_*.idat" -exec cp -t "$dest_dir" {} +
}

export -f copy_idats

# Use GNU Parallel to copy .idat files for each sample in parallel
cat "$sample_file" | xargs -P "$n_threads" -I {} bash -c "copy_idats '{}' '$destination_directory' '$raw_data_dir'"


echo "Copying valid sample IDATs completed."


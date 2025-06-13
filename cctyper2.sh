#!/bin/bash

# Check if exactly two arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

# Assign arguments to variables
input_dir="$1"
output_dir="$2"

# Check if input directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir" || {
    echo "Error: Failed to create output directory '$output_dir'."
    exit 1
}

# Check if input directory contains .fasta files
shopt -s nullglob # Prevent loop from running if no files match
fasta_files=("$input_dir"/*.fasta)
if [ ${#fasta_files[@]} -eq 0 ]; then
    echo "Error: No .fasta files found in '$input_dir'."
    exit 1
fi

# Activate conda environment
source /Users/pjoglekar/miniconda3/etc/profile.d/conda.sh
conda activate /Users/pjoglekar/miniforge3/envs/cctyper || {
    echo "Error: Failed to activate conda environment."
    exit 1
}

# Process each .fasta file
for f in "${fasta_files[@]}"; do
    base=$(basename "$f")
    out=$(echo "$base" | sed 's/_renamed.fasta/_cctyper_out/g')
    
    echo "Input file: $f"
    echo "Base name: $base"
    echo "Output name: $out"
    
    command="cctyper '$f' '$output_dir/$out' -t 20 --simplelog"
    echo "Command: $command"
    
    # Execute the command and append output to log
    cctyper "$f" "$output_dir/$out" -t 20 --simplelog >> "$output_dir/cctyper.log" 2>&1 || {
        echo "Error: cctyper failed for '$f'."
    }
done

# Deactivate conda environment
conda deactivate

echo "Processing complete. Output saved in '$output_dir'."
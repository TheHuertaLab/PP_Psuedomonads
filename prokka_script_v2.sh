#!/bin/bash
source /Users/pjoglekar/miniforge3/etc/profile.d/conda.sh # Load conda

conda activate prokka_env

# Define the directory containing your FASTA files
fasta_dir=$1
output_dir=$2

# Loop through each FASTA file in the directory
for fasta_file in "$fasta_dir"/*.fasta
do
  # Get the base name of the file (without the directory and extension)
  base_name=$(basename "$fasta_file" .fasta)
  
  # Run Prokka annotation
  prokka --outdir "$output_dir/$base_name" --prefix "$base_name" "$fasta_file" --cpus 15 --force
  mkdir -p "$output_dir"/prokka_gff_files

  for f in "$output_dir/$base_name"/*.gff; do
    mv "$f" "$output_dir"/prokka_gff_files
  done  
done

prokka --version
conda deactivate
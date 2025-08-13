#!/bin/bash

mkdir -p pseudomonas_renamed_genomes

dir=$1

for file in $dir/*.fasta; do

    echo "$file"
    base=$(basename "$file" | sed 's/.fasta//')
    echo "$base"

    out="$base" 
    echo "$out"
    
    awk '/^>/ {print ">'$base'_contig_" ++i; next} {print}' $file > pseudomonas_renamed_genomes/"$base".fasta

    
done
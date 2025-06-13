#!/bin/bash

mkdir -p phastest_renamed_contigs

dir=$1

for file in $dir/*.fasta; do

    echo "$file"
    base=$(basename "$file" | sed 's/.fasta//')
    echo "$base"

    out="$base"_renamed 
    echo "$out"
    
    awk '/^>/ {print ">'$base'_contig_" ++i; next} {print}' $file > phastest_renamed_contigs/"$base".fasta

    
done
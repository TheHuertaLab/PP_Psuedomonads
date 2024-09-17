#!/bin/bash

mkdir -p RenamedFastaFiles

dir=$1

for file in $dir/*.fasta; do

    echo "$file"
    base=$(basename "$file" | sed 's/.fasta//')

    out="$base"_renamed 
    echo "$out"
    
    awk '/^>/ {print ">'$base'_contig_" ++i; next} {print}' $file > RenamedFastaFiles/"$base"_renamed.fasta

    
done
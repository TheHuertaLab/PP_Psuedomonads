#!/bin/bash

source /Users/pjoglekar/miniconda3/etc/profile.d/conda.sh

conda activate /Users/pjoglekar/miniforge3/envs/cctyper

dir=$1

dir2=$2

mkdir -p $dir2

for f in $dir/*.fasta
do
  base=$(basename "$f")
  out=$(echo "$base" | sed 's/.fasta/_cctyper_out/g')
  
  echo "Input file: $f"
  echo "Base name: $base"
  echo "Output name: $out"
  
  command="cctyper $f $dir2/$out -t 10 --simplelog"
  echo "Command: $command"
  
  # Uncomment the line below if everything looks correct
  eval $command >> $dir2/cctyper.log 2>&1
done

conda deactivate

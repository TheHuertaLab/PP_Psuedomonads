#!/bin/bash

source /Users/pjoglekar/miniconda3/etc/profile.d/conda.sh

conda activate vibrant_env


dir=$1

for f in $dir/*.fasta
do
  base=$(basename "$f")
  out=$(echo "$base" | sed 's/.fasta/_out/g')
  
  echo "Input file: $f"
  echo "Base name: $base"
  echo "Output name: $out"
  
  command="python3 /Users/pjoglekar/miniconda3/envs/vibrant_env/bin/VIBRANT_run.py -t 10 -i $f -folder /Users/pjoglekar/work/pseudomonas/vibrant_out/$out"
  echo "Command: $command"
  
  # Uncomment the line below if everything looks correct
  eval $command 
done

conda deactivate

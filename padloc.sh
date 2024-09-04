#!/bin/bash

source /Users/pjoglekar/miniconda3/etc/profile.d/conda.sh

conda activate padloc
dir=$1

for file in $dir; 
do
mkdir $(basename $file | sed 's/.fasta/_padloc2_out/'); 
out=$(basename $file | sed 's/.fasta/_padloc2_out/');
padloc --fna $file --outdir $out --cpu 8;

done

conda deactivate

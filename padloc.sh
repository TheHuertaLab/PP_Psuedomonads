#!/bin/bash

source /Users/pjoglekar/miniforge3/etc/profile.d/conda.sh

conda activate /Users/pjoglekar/miniconda3/envs/padloc_env

dir=$1

for file in $dir/*.fna; 
do
mkdir $(basename $file | sed 's/.fna/_padloc2_out/'); 
out=$(basename $file | sed 's/.fna/_padloc2_out/');
padloc --fna $file --outdir $out --cpu 16;

done

conda deactivate

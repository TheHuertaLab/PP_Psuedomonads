#!/bin/bash

source /Users/pjoglekar/miniconda3/etc/profile.d/conda.sh

conda activate phispy_env

for file in /Users/pjoglekar/work/pseudomonas/pseudo_genbank/New_genomes/*.gb; 
do 
out=$(basename $file | sed 's/.gb/_out/');
phispy.py $file -o $out

done

conda deactivate

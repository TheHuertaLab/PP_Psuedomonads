#/bin/bash

#this scripts uses flye assembles to assemble pacbio reads in a given directory
#It can be run on a Daves_mac
#if using on another computer, you will have to chage path for conda.sh
#and change the path for the flye_env conda activate path/to/flye_env

#usage: bash flye.sh /path/to/directory/with/pacbio/reads. This will run in the terminal and,
##produce output in the terminal

#usage: nohup bash flye.sh /path/to/directory/with/pacbio/reads &. 
##This will run in the background. the std out and std error will be saved in nohup.out 
##in the directory where the command was executed

source /Users/pjoglekar/miniforge3/etc/profile.d/conda.sh # Load conda

conda activate flye_env #activates flye environment

dir=$1 

#path to fastq files from pacbio

echo "$dir"

mkdir -p $dir/flye_assemblies
chmod 777 $dir/flye_assemblies
out=$dir/flye_assemblies

for file in $dir/*.fastq.gz; do
    echo "$file"
    base=$(basename $file | sed 's/.fastq.gz//')

    flye --pacbio-hifi $file --out-dir $out/$base --genome-size 6m --plasmids  --threads 16 -i 5
done

flye --version
#!/bin/bash

mkdir -p renamed_phasted_files

for dir in /Users/pjoglekar/work/software/shared_docker_image/phastest-docker/phastest-app-docker/JOBS/*; do
    if [ -d "$dir" ]; then
        base=$(basename "$dir" | sed 's/_assembly//' | sed 's/_renamed/_phastest/')
        src_file="$dir/region_DNA.txt"
        dest_file="renamed_phasted_files/${base}.fasta"
        
        if [ -f "$src_file" ]; then
            cp "$src_file" "$dest_file"
            echo "Copied $src_file to $dest_file"
        else
            echo "File $src_file does not exist"
        fi
    fi
done

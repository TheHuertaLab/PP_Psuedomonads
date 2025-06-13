#!/bin/bash

cd /Users/pjoglekar/work/software/shared_docker_image/phastest-docker/phastest_inputs

for file in *.fasta; do
  FILE=$(basename "$file")
  echo "Processing file: $FILE"
  docker compose run phastest -i fasta --yes -s "$file" 
done
#!/bin/bash

cd /Users/pjoglekar/work/software/shared_docker_image/phastest-docker/phastest_inputs

for file in *.gbk; do
  # Extract the base name without the .gbk extension
  BASENAME=$(basename "$file" .gbk)
  echo "Processing file: $BASENAME.gbk"
  docker compose run phastest -i genbank -a "$BASENAME" --yes
done
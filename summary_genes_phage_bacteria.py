import os
import shutil

# Define the root directory
root_dir = "/Users/pjoglekar/work/pseudomonas/vibrant_out"

# Create the output directory
output_directory = "/Users/pjoglekar/work/pseudomonas/vib_summary_genomes"
os.makedirs(output_directory, exist_ok=True)

# Traverse the directory structure
for dirpath, dirnames, filenames in os.walk(root_dir):
    for filename in filenames:
        if filename.endswith("_combined.fna"):
            file_path = os.path.join(dirpath, filename)
            # Copy the file to the output directory
            shutil.copy(file_path, output_directory)

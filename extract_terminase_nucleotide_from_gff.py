from Bio import SeqIO
import os
import glob
import sys

# Get the directory path and output file name from command line arguments
directory = sys.argv[1]
output_file = sys.argv[2]

# Set to store unique names
unique_names = set()

# Iterate over all .gff files in the subdirectories
for file_path in glob.glob(os.path.join(directory, "**/*.gff"), recursive=True):
    # Open the file and extract the sequence for "terminase" and "large"
    with open(file_path, "r") as file, open(output_file, "a") as output:
        for record in SeqIO.parse(file, "gff"):
            if (
                "terminase" in record.description.lower()
                or "large" in record.description.lower()
            ):
                name = record.description
                if name not in unique_names:
                    unique_names.add(name)
                    # Extract the sequence from the appropriate feature in the GFF record
                    sequence = record.seq
                    output.write(f">{name}\n{sequence}\n")

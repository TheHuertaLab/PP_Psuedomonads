import os
import sys
import pandas as pd

# Directory path where the .fasta files are located
directory = sys.argv[1]

# Create a new file to write the output
output_file = sys.argv[2]
with open(output_file, "w") as f_out:

    # Iterate over all files in the directory
    f_out.write("accession_no\tcontig_id\n")
    for filename in os.listdir(directory):
        if filename.endswith(".fasta"):
            # Get the genome name (the first part of the file name)
            genome_name = filename.split(".")[0]

            # Open the file and read the contig names
            contig_names = set()  # use a set to store unique contig names
            with open(os.path.join(directory, filename), "r") as file:
                for line in file:
                    if line.startswith(">"):
                        # extract the general contig id
                        contig_id = line.strip()[1:].split("_")[0]
                        contig_names.add(contig_id)

            # Write the genome name and contig names to the output file
            f_out.write(f"{genome_name}\t{':'.join(contig_names)}\n")

# Print a message indicating the output file path
print(f"Output file created: {output_file}")

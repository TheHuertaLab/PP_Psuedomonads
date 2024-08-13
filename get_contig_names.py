import os
import sys
import pandas as pd

# Directory path where the .fasta files are located
directory = sys.argv[1]

# Create a new file to write the output
output_file = sys.argv[2]

# second directory with pseudomonas type strain sequences
# dir2 = "/Users/pjoglekar/work/pseudomonas/pseudo_fasta"

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


##### This code will add contig names from the reference sequences to output file
# with open(output_file, "a") as f_out:
#     serial_number = 1  # Start the serial number from 1

#     for filename in os.listdir(dir2):
#         if filename.endswith(".fasta"):
#             # Replace genome name with serial number
#             genome_name = str(serial_number)

#             # Open the file and read the contig names
#             contig_names = set()  # Use a set to store unique contig names
#             with open(os.path.join(dir2, filename), "r") as file:
#                 for line in file:
#                     if line.startswith(">"):
#                         # Extract the general contig id
#                         contig_id = line.strip()[1:].split(" ")[0]
#                         contig_names.add(contig_id)

#             # Write the serial number and contig names to the output file, each on a new line
#             for contig in contig_names:
#                 f_out.write(f"{genome_name}\t{contig}\n")

#             # Increment the serial number for the next file
#             serial_number += 1

# Print a message indicating the output file path
print(f"Output file created: {output_file}")

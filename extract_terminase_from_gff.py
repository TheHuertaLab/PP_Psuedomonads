from Bio import SeqIO
import os

directory = "/Users/pjoglekar/work/Repre_Genomes/vib_summary_faa"  # Replace with the actual directory path

# Create a new file to write the sequences
output_file = os.path.join(directory, "all_terminase7.faa")

# Set to store unique names
unique_names = set()

# Iterate over all files in the directory
for filename in os.listdir(directory):
    if filename.endswith("combined.faa"):
        file_path = os.path.join(directory, filename)

        # Open the file and extract the sequence for "terminase"
        with open(file_path, "r") as file, open(output_file, "a") as output:
            for record in SeqIO.parse(file, "fasta"):
                if (
                    "terminase" in record.description
                    or "Terminase" in record.description
                ) and "small" not in record.description:
                    name = record.description
                    if name not in unique_names:
                        unique_names.add(name)
                        output.write(f">{name}\n{record.seq}\n")

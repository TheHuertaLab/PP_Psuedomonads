from Bio import SeqIO
import os

directory = "/Users/pjoglekar/work/pseudomonas/pseudo_genomes/pseudo_cucurbit_genomes"  # Replace with the actual directory path

# Create a new file to write the sequences
output_file = os.path.join(directory, "all_gyrB.fna")

# Iterate over all files in the directory and its subdirectories
for root, dirs, files in os.walk(directory):
    for filename in files:
        if filename.endswith("genes.fna"):
            file_path = os.path.join(root, filename)

            # Open the file and extract the sequences for "DNA gyrase B"
            with open(file_path, "r") as file, open(output_file, "a") as output:
                for record in SeqIO.parse(file, "fasta"):
                    with open(
                        "/Users/pjoglekar/work/pseudomonas/pseudo_genomes/list_nucleotides_location_gyrB.txt",
                        "r",
                    ) as desc_file:
                        for line in desc_file:
                            if line.strip() in record.description:
                                name = record.description
                                output.write(f">{name}\n{record.seq}\n")

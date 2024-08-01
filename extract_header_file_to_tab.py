import os
from Bio import SeqIO

input_file = "/Users/pjoglekar/work/Repre_Genomes/vib_summary_faa/all_terminase7.faa"
output_file = "/Users/pjoglekar/work/Repre_Genomes/vib_summary_faa/all_terminase7.tsv"

with open(input_file, "r") as file, open(output_file, "w") as output:
    for record in SeqIO.parse(file, "fasta"):
        record_description = record.description
        tab_separated = "\t".join([record_description])
        output.write(tab_separated + "\n")

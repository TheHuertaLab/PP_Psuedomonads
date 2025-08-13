import os
from glob import glob
from Bio import SeqIO

def get_protein_info_from_prodigal_gff(prodigal_gff, start, end):
    with open(prodigal_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            gff_start, gff_end = int(fields[3]), int(fields[4])
            if gff_start == start and gff_end == end:
                # Extract protein number from attributes
                protein_number = None
                for attr in fields[8].split(";"):
                    if attr.startswith("ID="):
                        protein_number = attr.split("_")[-1]
                return protein_number
    return None

def extract_matching_proteins(subfolder):
    padloc_gff = glob(os.path.join(subfolder, "*padloc.gff"))[0]
    prodigal_gff = glob(os.path.join(subfolder, "*prodigal.gff"))[0]
    prodigal_faa = glob(os.path.join(subfolder, "*prodigal.faa"))[0]
    
    # Add directory basename as prefix to output file
    basename = os.path.basename(subfolder)
    output_faa = os.path.join(subfolder, f"{basename}_padloc_matched_proteins.faa")

    # Step 1: Get all (start, end) and Name from padloc.gff
    matches = []
    with open(padloc_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            start, end = int(fields[3]), int(fields[4])
            name = None
            for attr in fields[8].split(";"):
                if attr.startswith("Name="):
                    name = attr.split("=")[-1]
            matches.append((start, end, name))

    # Step 2: For each (start, end), get protein number from prodigal.gff
    protein_info = []
    for start, end, name in matches:
        num = get_protein_info_from_prodigal_gff(prodigal_gff, start, end)
        if num:
            protein_info.append((num, name))

    # Step 3: Extract proteins from prodigal.faa and add padloc name to FASTA header
    with open(output_faa, "w") as out_faa:
        for record in SeqIO.parse(prodigal_faa, "fasta"):
            for num, name in protein_info:
                if record.id.endswith(f"_{num}"):
                    # Modify the FASTA header to include the padloc name
                    record.id = f"{record.id} | padloc_name={name}" if name else record.id
                    record.description = ""  # Clear description
                    SeqIO.write(record, out_faa, "fasta")

# Main: process all subfolders
parent_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new"  # Change to your parent directory
for subfolder in glob(os.path.join(parent_dir, "*")):
    if os.path.isdir(subfolder):
        try:
            extract_matching_proteins(subfolder)
            print(f"Processed {subfolder}")
        except Exception as e:
            print(f"Error in {subfolder}: {e}")
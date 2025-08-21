import os
from glob import glob
from Bio import SeqIO

def get_protein_number_from_prodigal_gff(prodigal_gff, name, start, end):
    with open(prodigal_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            gff_name, gff_start, gff_end = int(fields[0]), int(fields[3]), int(fields[4])
            if gff_name == name, gff_start == start and gff_end == end:
                # Extract protein number from attributes (usually ID=..._number)
                for attr in fields[8].split(";"):
                    if attr.startswith("ID="):
                        return attr.split("_")[-1]
    return None

def extract_matching_proteins(subfolder):
    padloc_gff = glob(os.path.join(subfolder, "*padloc.gff"))[0]
    prodigal_gff = glob(os.path.join(subfolder, "*prodigal.gff"))[0]
    prodigal_faa = glob(os.path.join(subfolder, "*prodigal.faa"))[0]
    output_faa = os.path.join(subfolder, "padloc_matched_proteins.faa")

    # Step 1: Get all (start, end) from padloc.gff
    matches = []
    with open(padloc_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            name, start, end = int(fields[0]),int(fields[3]), int(fields[4])
            matches.append((name, start, end))

    # Step 2: For each (start, end), get protein number from prodigal.gff
    protein_numbers = []
    for name, start, end in matches:
        num = get_protein_number_from_prodigal_gff(prodigal_gff, name, start, end)
        if num:
            protein_numbers.append(num)

    # Step 3: Extract proteins from prodigal.faa
    with open(output_faa, "w") as out_faa:
    for record in SeqIO.parse(prodigal_faa, "fasta"):
        # Extract contig and protein number from record.id
        # Example: Ps_papulans_renamed_contig_35_39
        parts = record.id.split("_")
        record_contig = "_".join(parts[:-1])
        record_num = parts[-1]
        for contig, num, name in protein_info:
            if record_contig == contig and record_num == num:
                record.id = f"{record.id} | padloc_name={name}" if name else record.id
                record.description = ""
                SeqIO.write(record, out_faa, "fasta")

# Main: process all subfolders
parent_dir = "/path/to/parent"  # Change to your parent directory
for subfolder in glob(os.path.join(parent_dir, "*")):
    if os.path.isdir(subfolder):
        try:
            extract_matching_proteins(subfolder)
            print(f"Processed {subfolder}")
        except Exception as e:
            print(f"Error in {subfolder}: {e}")
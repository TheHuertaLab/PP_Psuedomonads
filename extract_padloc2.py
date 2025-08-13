import os
from glob import glob
from Bio import SeqIO

def get_protein_number_from_prodigal_gff(prodigal_gff, contig, start, end):
    with open(prodigal_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            gff_contig = fields[0]
            gff_start, gff_end = int(fields[3]), int(fields[4])
            gff_id = None
            for attr in fields[8].split(";"):
                if attr.startswith("ID="):
                    gff_id = attr.split("=")[-1]
            if gff_contig == contig and gff_start == start and gff_end == end and gff_id:
                protein_num = gff_id.split("_")[-1]
                return protein_num
    return None

def extract_matching_proteins(subfolder):
    padloc_gff = glob(os.path.join(subfolder, "*padloc.gff"))[0]
    prodigal_gff = glob(os.path.join(subfolder, "*prodigal.gff"))[0]
    prodigal_faa = glob(os.path.join(subfolder, "*prodigal.faa"))[0]
    
    basename = os.path.basename(subfolder)
    output_faa = os.path.join(subfolder, f"{basename}_padloc_matched_proteins.faa")

    # Step 1: Get all (contig, start, end, name) from padloc.gff
    matches = []
    with open(padloc_gff) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            contig = fields[0]
            start, end = int(fields[3]), int(fields[4])
            name = None
            for attr in fields[8].split(";"):
                if attr.startswith("Name="):
                    name = attr.split("=")[-1]
            matches.append((contig, start, end, name))

    # Step 2: For each (contig, start, end, name), get protein number from prodigal.gff
    protein_info = []
    for contig, start, end, name in matches:
        num = get_protein_number_from_prodigal_gff(prodigal_gff, contig, start, end)
        if num:
            protein_info.append((contig, num, name))

    # Step 3: Extract proteins from prodigal.faa and add padloc name to FASTA header
    with open(output_faa, "w") as out_faa:
        for record in SeqIO.parse(prodigal_faa, "fasta"):
            idx = record.id.rfind("_")
            record_contig = record.id[:idx]
            record_num = record.id[idx+1:]
            for contig, num, name in protein_info:
                if record_contig == contig and record_num == num:
                    record.id = f"{record.id} | padloc_name={name}" if name else record.id
                    record.description = ""
                    SeqIO.write(record, out_faa, "fasta")
                    break  # Only write once per match

# Main: process all subfolders
parent_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new"
for subfolder in glob(os.path.join(parent_dir, "*")):
    if os.path.isdir(subfolder):
        try:
            extract_matching_proteins(subfolder)
            print(f"Processed {subfolder}")
        except Exception as e:
            print(f"Error processing {subfolder}: {e}")
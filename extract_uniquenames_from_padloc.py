import os
from glob import glob

base_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new"
output_file = os.path.join(base_dir, "unique_padloc_names.txt")

unique_names = set()

# Find all *_padloc_matched_proteins.faa files in all subdirectories
faa_files = glob(os.path.join(base_dir, "*", "*_padloc_matched_proteins.faa"))

for faa_file in faa_files:
    with open(faa_file) as f:
        for line in f:
            if line.startswith(">"):
                # Find the padloc_name in header
                parts = line.strip().split("padloc_name=")
                if len(parts) > 1:
                    name = parts[1].strip()
                    unique_names.add(name)

# Save unique names to file
with open(output_file, "w") as out:
    for name in sorted(unique_names):
        out.write(name + "\n")

print(f"✅ Extracted {len(unique_names)} unique padloc names to {output_file}")

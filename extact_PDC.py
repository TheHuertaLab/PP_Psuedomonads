import os

def extract_pdc_proteins_from_all_files(parent_dir, output_dir, pdc_protein_names):
    """
    Extract sequences for each PDC protein from all matching .faa files in the parent directory and its subdirectories.

    :param parent_dir: Parent directory containing subdirectories with *_padloc_matched_proteins.faa files.
    :param output_dir: Directory to save the output .faa files.
    :param pdc_protein_names: List of PDC protein names to extract.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Initialize a dictionary to store sequences for each PDC protein
    pdc_sequences = {pdc: [] for pdc in pdc_protein_names}

    # Walk through all subdirectories in the parent directory
    for root, _, files in os.walk(parent_dir):
        for file in files:
            if file.endswith("_padloc_matched_proteins.faa"):
                input_file = os.path.join(root, file)
                with open(input_file, 'r') as f:
                    current_header = None
                    current_sequence = []
                    
                    for line in f:
                        line = line.strip()
                        if line.startswith(">"):  # Header line
                            # Save the previous sequence if it matches a PDC protein
                            if current_header and current_sequence:
                                for pdc in pdc_protein_names:
                                    if f"padloc_name={pdc}" in current_header:
                                        pdc_sequences[pdc].append(current_header)
                                        pdc_sequences[pdc].append("".join(current_sequence))
                                        break
                            
                            # Start a new sequence
                            current_header = line
                            current_sequence = []
                        else:
                            # Sequence line
                            current_sequence.append(line)
                    
                    # Save the last sequence
                    if current_header and current_sequence:
                        for pdc in pdc_protein_names:
                            if f"padloc_name={pdc}" in current_header:
                                pdc_sequences[pdc].append(current_header)
                                pdc_sequences[pdc].append("".join(current_sequence))
                                break

    # Write sequences to separate files
    for pdc, sequences in pdc_sequences.items():
        if sequences:  # Only create files for PDC proteins with sequences
            output_file = os.path.join(output_dir, f"{pdc}.faa")
            with open(output_file, 'w') as out_file:
                out_file.write("\n".join(sequences))
                out_file.write("\n")

# Define the parent directory and output directory
parent_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new"
output_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new/PDC_proteins"

# List of PDC protein names
pdc_protein_names = [
    "PDC-S04", "PDC-S15", "PDC-S21", "PDC-M01A", "PDC-M01B","PDC-S65", "PDC-S35", "PDC-S08", "PDC-S06", "PDC-S12", "PDC-S28",
    "PDC-S39", "PDC-S02", "PDC-M32A","PDC-M32B", "PDC-S24", "PDC-S64", "PDC-M54A","PDC-M54B", "PDC-S37", "PDC-S45", "PDC-S32", "PDC-S47",
    "PDC-S52", "PDC-S58", "PDC-S44", "PDC-M30A","PDC-M30B", "PDC-M03A","PDC-M03A","PDC-S31", "PDC-S25", "PDC-M55", "PDC-M65A","PDC-M65B", "PDC-M66A","PDC-M66B",
    "PDC-M02", "PDC-M38A","PDC-M38B", "PDC-S73", "PDC-S33", "PDC-S29", "PDC-M18A", "PDC-M18B","PDC-S05", "PDC-S36", "PDC-S71", "PDC-M24A","PDC-M24B",
    "PDC-M07A", "PDC-M07B","PDC-S51", "PDC-S16", "PDC-S63"
]

# Run the extraction
extract_pdc_proteins_from_all_files(parent_dir, output_dir, pdc_protein_names)
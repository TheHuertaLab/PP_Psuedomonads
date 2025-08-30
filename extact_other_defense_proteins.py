import os

def extract_proteins_from_all_files(parent_dir, output_dir, defense_protein_names):
    """
    Extract sequences for each defense protein from all matching .faa files in the parent directory and its subdirectories.

    :param parent_dir: Parent directory containing subdirectories with *_padloc_matched_proteins.faa files.
    :param output_dir: Directory to save the output .faa files.
    :param defense_protein_names: List of defense protein names to extract.
    """
    os.makedirs(output_dir, exist_ok=True)
    defense_sequences = {defense: [] for defense in defense_protein_names}

    for root, _, files in os.walk(parent_dir):
        for file in files:
            if file.endswith("_padloc_matched_proteins.faa"):
                input_file = os.path.join(root, file)
                with open(input_file, 'r') as f:
                    current_header = None
                    current_sequence = []
                    for line in f:
                        line = line.rstrip()
                        if line.startswith(">"):
                            if current_header and current_sequence:
                                for defense in defense_protein_names:
                                    if f"padloc_name={defense}" in current_header:
                                        defense_sequences[defense].append(current_header + "\n" + "".join(current_sequence))
                                        break
                            current_header = line
                            current_sequence = []
                        else:
                            current_sequence.append(line)
                    if current_header and current_sequence:
                        for defense in defense_protein_names:
                            if f"padloc_name={defense}" in current_header:
                                defense_sequences[defense].append(current_header + "\n" + "".join(current_sequence))
                                break

    for defense, sequences in defense_sequences.items():
        if sequences:
            output_file = os.path.join(output_dir, f"{defense}.faa")
            with open(output_file, 'w') as out_file:
                out_file.write("\n".join(sequences) + "\n")



#define parent and output directories
parent_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new"
output_dir = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/pacbio_genomes/padloc_new/defense_proteins"

defense_protein_names = [
    "ATPase-Toprim_I-B", "AbiD", "AbiEi", "AbiEii", "AbiLi", "AbiLii", "AbiO",
    "AbiU", "Aga_SIR2-PNTD1-PNTD2", "AriA", "AriB", "Avs2", "Avs4", "BrxA",
    "BrxB", "BrxC", "BrxHI", "BrxL", "CSD_V", "Cas1f", "Cas23f",
    "Cas5f", "Cas6f", "Cas7f", "Cas8f", "Control_protein", "Cyclase", "DUF1156_MTase",
    "DUF2290", "DUF3780", "DUF4238-Pers", "DUF4297", "DUF499", "DarG", "DarT",
    "DndA", "DndB", "DndC", "DndD", "DndE", "DpdA", "DpdB",
    "DpdC", "DpdD", "DpdE", "DpdF", "DpdG", "DpdH", "DpdJ",
    "DpdK", "DrmA", "DrmB", "DrmC", "DrmD", "DrmMI", "DrmMII",
    "Drt1a", "Drt2", "Drt3a", "Drt3b", "Drt4", "DruA1", "DruB1",
    "DruC1", "DruD1", "DruE1", "Dsr1", "DzbA", "DzbB", "E1-E2",
    "Effector", "GajA", "GajB", "HEC-01A", "HEC-01B", "HEC-02A", "HEC-02B",
    "HEC-06", "HEC-09", "HORMA", "HTH_VI", "HamA1", "HamB1", "Helicase",
    "HerA", "Hydrolase", "Hydrolase-TM", "IetA", "IetS", "JAB", "JetA1",
    "JetB1", "JetC1", "JetD1", "KwaA", "KwaB", "LeoA", "LeoBC",
    "LmuA", "LmuB", "LmuC", "MTase_I", "MTase_II", "MTase_III", "MkoA",
    "MkoB", "MkoC", "MzaB", "MzaC", "MzaD", "MzaE", "NDT_II-A",
    "Nhi", "NsnA", "NsnB", "NsnC", "Old", "Orf5", "PD-Lambda-5_A",
    "PD-Lambda-5_B", "PD-T4-5", "PD-T4-6", "PD-T7-1", "PD-T7-3", "PD-T7-5", "PLD_Helicase",
    "PP2C", "PRTase_III-A", "PbeA", "PglX", "PglZ", "PprDinG", "PprPIWI",
    "PprRE", "PrrC", "PsyrA", "PsyrT", "PtuA1", "PtuB1", "PycC",
    "PycTM", "QatA", "QatB", "QatC", "QatD", "REase_I",
    "REase_II", "REase_III", "REase_MTase_IIG", "RT-Protease_XI", "RT-Toprim_I-C", "RT_I-B", "RT_II-A",
    "RT_III-A", "RT_UG17", "RT_UG5-nitrilase", "RT_UG9_small", "RT_V", "RT_VI", "RTasNERD-Hel",
    "RTasPolA", "RTasSLATT", "RmrA", "RmrT", "SIR2", "STK_OB", "SduA",
    "ShosA", "ShosT", "SngA", "SngB", "SngC", "SoFic", "Specificity_I",
    "SspB", "SspC", "SspD", "SspF", "SspG", "SspH", "TIR-NLR",
    "TRIP13", "TerY", "ThsA1", "ThsB1", "Tmn", "Upx", "UzuA",
    "Vsr", "WYL", "ZacA", "ZacB", "ZacC", "bSEFIR", "capH",
    "capP", "dGTPase", "mREase_II", "mREase_IV", "mad1", "mad3", "mad5",
    "mad8", "pAgo_I"
]

extract_proteins_from_all_files(parent_dir, output_dir, defense_protein_names)


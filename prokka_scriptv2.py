#!/Users/pjoglekar/miniforge3/bin/python3

import os
import subprocess
import sys
from glob import glob
from pathlib import Path
import shutil

def main(fasta_dir, output_dir):
    # Prefix for conda activation
    conda_prefix = "source /Users/pjoglekar/miniforge3/etc/profile.d/conda.sh && conda activate prokka_env && "

    fasta_files = glob(os.path.join(fasta_dir, "*.fasta"))
    gff_dir = os.path.join(output_dir, "prokka_gff_files")
    os.makedirs(gff_dir, exist_ok=True)

    for fasta_file in fasta_files:
        base_name = Path(fasta_file).stem
        out_path = os.path.join(output_dir, base_name)
        os.makedirs(out_path, exist_ok=True)
        # Run Prokka with conda activation
        prokka_cmd = (
            f"{conda_prefix}prokka --outdir '{out_path}' --prefix '{base_name}' '{fasta_file}' --cpus 15 --force"
        )
        subprocess.run(prokka_cmd, shell=True, executable="/bin/bash", check=True)
        # Move GFF files with conda activation (mv)
        for gff_file in glob(os.path.join(out_path, "*.gff")):
            mv_cmd = f"{conda_prefix}mv '{gff_file}' '{gff_dir}'"
            subprocess.run(mv_cmd, shell=True, executable="/bin/bash", check=True)

    # Print Prokka version with conda activation
    subprocess.run(f"{conda_prefix}prokka --version", shell=True, executable="/bin/bash")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python prokka_script_v2.py <fasta_dir> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
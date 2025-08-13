#!/Users/pjoglekar/miniforge3/bin/python3

import os
import subprocess
import sys
from glob import glob
from pathlib import Path
import shutil

def main(fasta_dir,output_dir):
    # Prefix for conda activation
    conda_prefix = "source /Users/pjoglekar/miniforge3/etc/profile.d/conda.sh && conda activate /Users/pjoglekar/miniconda3/envs/padloc_env &&"
    
    fasta_files = glob(os.path.join(fasta_dir, "*.fasta"))
    
    for fasta_file in fasta_files:
        base_name = Path(fasta_file).stem
        out_path = os.path.join(output_dir, base_name)
        os.makedirs(out_path, exist_ok=True)
        #Run padloc command with conda activation
        padloc_cmd = (
            f"{conda_prefix}padloc --fna {fasta_file} --outdir {out_path} --cpu 16"
        )
        subprocess.run(padloc_cmd, shell=True, executable="/bin/bash", check=True)
        
    #Print padloc version
    subprocess.run(f"{conda_prefix}padloc --version", shell=True, executable="/bin/bash", check=True)
    
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python padloc.py <fasta_dir> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
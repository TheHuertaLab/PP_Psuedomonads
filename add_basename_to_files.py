import os

#flye produces files with assembly.fasta name for each assembly and so on. 
# This script adds the basename of the directory to the file name.
def add_basename_to_files(root_dir):
    for subdir, _, files in os.walk(root_dir):
        basename = os.path.basename(subdir)
        for file_name in files:
            old_file_path = os.path.join(subdir, file_name)
            new_file_name = f"{basename}_{file_name}"
            new_file_path = os.path.join(subdir, new_file_name)
            os.rename(old_file_path, new_file_path)
            print(f"Renamed: {old_file_path} to {new_file_path}")

if __name__ == "__main__":
    root_directory = "/Users/pjoglekar/work/reads/Feb_2025/SQ031_Huerta_pool3_04152025/fastq_files/flye_assemblies"
    add_basename_to_files(root_directory)
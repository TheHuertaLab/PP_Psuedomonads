import os
import shutil

#sometimes you have to move things around. This script will help you copy files 
##from one directory to another based on a list of number of names you are looking for

#change input and output folder accordingly

input_folder = "/Users/pjoglekar/work/pseudomonas/pseudo_prokka_out/prokka_gff_files"
output_folder = "/Users/pjoglekar/work/pseudomonas/Pseudomonas_data_from_selected_genomes/gff_files"

#Numbers to looks for in the directory 
#I used numbers, but it could be any string

numbers = [
    "2990935028", "2990761813", "2990977861", "2991173094", "2991383323", 
    "2990827573", "2991049539", "2993069567", "2992283825", "2992154122", 
    "2992789122", "2991837845", "2992832258", "2992262078", "2992555032", 
    "2992143505", "2992638703", "2992007228", "2992633404", "2992837609", 
    "2992816361", "2992098980", "2991765579", "2991365836", "2992924922", 
    "2992605350"
]

#check if output folder exists, if not create it
os.makedirs(output_folder,exist_ok=True)

#Iterate through the folder and copy files
for file in os.listdir(input_folder):
    #Check if the filei s a .fasta file and starts with one of the numbers in the above list
    if file.endswith(".gff") and any(file.startswith(num) for num in numbers):
        #Copy the file to the output folder
        shutil.copy(os.path.join(input_folder,file),os.path.join(output_folder,file))
        print(f"Copying {file} to {output_folder}")
        
print("Done")
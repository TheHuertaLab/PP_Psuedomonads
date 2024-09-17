import os

# Define the input and output paths
jobs_path = "/Users/pjoglekar/work/software/shared_docker_image/phastest-docker/phastest-app-docker/JOBS/"
output_file = "/Users/pjoglekar/work/software/shared_docker_image/phastest-docker/phastest-app-docker/JOBS/region_DNA.txt"

# Function to process the sequence files in each directory
def process_txt_files(dir_path, region_name):
    for file_name in os.listdir(dir_path):
        file_path = os.path.join(dir_path, file_name)
        
        # Only process .txt files
        if os.path.isfile(file_path) and file_name.endswith(".txt"):
            # Open the .txt file and modify the headers
            with open(file_path, 'r') as txt_file:
                lines = txt_file.readlines()
            
            # Rewrite the file with modified headers
            with open(file_path, 'w') as txt_file:
                for line in lines:
                    if line.startswith('>'):
                        # Check if the region_name appears multiple times
                        if line.startswith(f">{region_name}_"):
                            # Remove extra occurrences of region_name
                            new_header = line.replace(f"{region_name}_", '', line.count(f"{region_name}_") - 1)
                            txt_file.write(new_header)
                        else:
                            # Add region_name only if not already present
                            new_header = f">{region_name}_{line[1:]}"
                            txt_file.write(new_header)
                    else:
                        txt_file.write(line)

# List to track region names we've already written to avoid duplicates
written_regions = []

# Read the region_DNA.txt file if it already exists to avoid overwriting existing content
if os.path.exists(output_file):
    with open(output_file, 'r') as outfile:
        written_regions = [line.strip() for line in outfile.readlines()]

# Open the output file for writing the region names
with open(output_file, 'a') as outfile:
    # Iterate through each subdirectory in the JOBS path
    for dir_name in os.listdir(jobs_path):
        dir_path = os.path.join(jobs_path, dir_name)
        
        # Check if it is a directory
        if os.path.isdir(dir_path):
            # Split and extract the part before the first '_'
            region_name = dir_name.split('_')[0]

            # If this region_name hasn't already been written, write it
            if region_name not in written_regions:
                outfile.write(region_name + '\n')
                written_regions.append(region_name)  # Add to the list to track

            # Process the .txt files in this directory
            process_txt_files(dir_path, region_name)

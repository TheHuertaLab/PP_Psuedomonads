import os
import sys
import pandas as pd

# Get the directory path and output file path from command line arguments
directory_path = sys.argv[1]
output_file_path = sys.argv[2]

# Initialize an empty DataFrame to store the concatenated data
concatenated_data = pd.DataFrame()

# Iterate over all subfolders in the directory
for root, dirs, files in os.walk(directory_path):
    for file in files:
        if file.endswith(".csv"):
            # Read the CSV file
            file_path = os.path.join(root, file)
            data = pd.read_csv(file_path)

            # Keep only the desired columns
            data = data[["seqid", "system"]]

            # Concatenate the data to the main DataFrame
            concatenated_data = pd.concat([concatenated_data, data])

# Write the concatenated data to the output file
concatenated_data.to_csv(output_file_path, index=False)

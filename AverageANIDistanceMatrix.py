import pandas as pd
import numpy as np

# Load the CSV file into a pandas DataFrame
file_path = '/Users/pjoglekar/work/prem/ralstonia_genomes_letters/ani_out/ani_out_new/ANIm_percentage_identity.tab'  # Replace with your actual CSV file path
ani_df = pd.read_csv(file_path, sep='\t', index_col=0)

# Get the matrix values as a NumPy array
ani_matrix = ani_df.values

# Get the number of rows/columns (since it's a square matrix)
n = ani_matrix.shape[0]

# Create a new matrix to store the averaged values
averaged_matrix = np.zeros_like(ani_matrix)

# Loop through the matrix and calculate the average for each pairwise comparison
for i in range(n):
    for j in range(i, n):  # Only iterate through upper triangular matrix including diagonal
        avg_value = (ani_matrix[i, j] + ani_matrix[j, i]) / 2
        averaged_matrix[i, j] = avg_value
        averaged_matrix[j, i] = avg_value

# Create a new DataFrame for the averaged matrix, preserving the original row/column labels
averaged_ani_df = pd.DataFrame(averaged_matrix, index=ani_df.index, columns=ani_df.columns)

# Save the averaged matrix to a new CSV file
output_file_path = '/Users/pjoglekar/work/prem/ralstonia_genomes_letters/ani_out/ani_out_new/corrected_ANIm_percentage_identity.tab'
averaged_ani_df.to_csv(output_file_path, sep='\t')

# Return the output file path for downloading
output_file_path
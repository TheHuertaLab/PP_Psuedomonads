import os

# Get the list of files in the directory
files = os.listdir("/Users/pjoglekar/work/pseudomonas/vib_summary_genomes")

# Iterate over the files
for filename in files:
    # Check if the filename starts with '2'
    if filename.startswith("2"):
        # Add 'Pseudomonas_' to the filename
        new_filename = "Pseudomonas_" + filename

        # Rename the file
        os.rename(
            os.path.join(
                "/Users/pjoglekar/work/pseudomonas/vib_summary_genomes", filename
            ),
            os.path.join(
                "/Users/pjoglekar/work/pseudomonas/vib_summary_genomes", new_filename
            ),
        )

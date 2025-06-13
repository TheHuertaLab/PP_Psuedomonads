import os
import csv
import gzip

# Path to the CSV file containing barcode-to-sample mapping
csv_file = "/Users/pjoglekar/work/reads/Feb_2025/barcodes.csv"
fastq_dir = "/Users/pjoglekar/work/reads/Feb_2025/SQ031_Huerta_pool3_04152025/fastq_files"

# Verify paths exist
if not os.path.exists(csv_file):
    print(f"Error: CSV file {csv_file} does not exist")
    exit(1)
if not os.path.exists(fastq_dir):
    print(f"Error: Directory {fastq_dir} does not exist")
    exit(1)

# Read the barcode-to-sample mapping
barcode_to_sample = {}
with open(csv_file, mode='r') as file:
    reader = csv.DictReader(file)
    if 'barcode' not in reader.fieldnames or 'sample_name' not in reader.fieldnames:
        print("Error: CSV must have 'barcode' and 'sample_name' columns")
        exit(1)
    for row in reader:
        barcode_to_sample[row['barcode']] = row['sample_name']
print(f"Loaded {len(barcode_to_sample)} barcode-to-sample mappings")

# Rename the FASTQ files
found_files = False
for filename in os.listdir(fastq_dir):
    if filename.endswith(".fastq.gz"):
        found_files = True
        parts = filename.split('.')
        print(f"Processing file: {filename}, parts: {parts}")
        if len(parts) > 2:
            barcode_full = parts[2] #here it takes the full regions bcXXXX--bcXXXX
            barcode = barcode_full.split("--")[0] 
            if barcode in barcode_to_sample:
                sample_name = barcode_to_sample[barcode]
                new_filename = f"{sample_name}.fastq.gz"
                old_path = os.path.join(fastq_dir, filename)
                new_path = os.path.join(fastq_dir, new_filename)
                try:
                    os.rename(old_path, new_path)
                    print(f"Renamed: {filename} -> {new_filename}")
                except OSError as e:
                    print(f"Error renaming {filename}: {e}")
            else:
                print(f"Barcode {barcode} not found in CSV")
        else:
            print(f"Filename {filename} does not have expected format")
if not found_files:
    print(f"No .fastq.gz files found in {fastq_dir}")
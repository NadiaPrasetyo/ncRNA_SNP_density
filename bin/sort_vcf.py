import pandas as pd

# Function to check if the chromosome is between 1 and 22
def is_valid_chromosome(chromosome):
    # Extract the numeric part of the chromosome name (e.g., "chr1", "chr2", ...)
    chrom_name = chromosome[3:]  # Skip the "chr" prefix
    if chrom_name.isdigit():
        # Include only chromosomes 1-22
        return 1 <= int(chrom_name) <= 22
    # Exclude non-numeric chromosomes (e.g., "chrX", "chrY", "chrM")
    return False

# Function to handle chromosome sorting
def parse_chromosome(chromosome):
    # Try to convert the chromosome number to an integer, otherwise return the string
    if chromosome[3:].isdigit():
        return int(chromosome[3:])
    return chromosome

# Define the path to your VCF file
input_file = 'data/decomposed.vcf'
output_file = 'data/sorted_all_variants.vcf'

# Read the VCF file, skipping the header lines
with open(input_file, 'r') as infile:
    lines = infile.readlines()

# Separate the header and data lines
header_lines = []
data_lines = []

for line in lines:
    if line.startswith('#'):
        header_lines.append(line)
    else:
        data_lines.append(line)

# Filter out lines where REF or ALT are more than 40 characters
# Also, exclude chromosomes beyond 1-22
filtered_data_lines = [
    line for line in data_lines
    if len(line.split('\t')[3]) <= 40 and len(line.split('\t')[4]) <= 40
    and is_valid_chromosome(line.split('\t')[0])
]

# Sort the filtered data lines by chromosome number
sorted_data_lines = sorted(filtered_data_lines, key=lambda x: parse_chromosome(x.split('\t')[0]))

# Write the sorted and filtered VCF file back with the original header
with open(output_file, 'w') as outfile:
    outfile.writelines(header_lines)
    outfile.writelines(sorted_data_lines)

print(f'Sorted and filtered VCF file has been saved to {output_file}')

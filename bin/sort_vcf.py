import pandas as pd

# Define the path to your VCF file
input_file = 'data/filtered_variants.vcf'
output_file = 'data/sorted_filtered_variants.vcf'

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

# Sort the data lines by chromosome number
# Since the chromosome number is the first part of the line (before the first tab)
sorted_data_lines = sorted(data_lines, key=lambda x: (int(x.split('\t')[0][3:]) if x.split('\t')[0][3:].isdigit() else x.split('\t')[0]))

# Write the sorted VCF file back with the original header
with open(output_file, 'w') as outfile:
    outfile.writelines(header_lines)
    outfile.writelines(sorted_data_lines)

print(f'Sorted VCF file has been saved to {output_file}')

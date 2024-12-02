import re
import csv
from collections import defaultdict, Counter

# Function to process the text file and extract data
def process_file(input_file, output_file):
    gene_data = defaultdict(Counter)
    
    # Open and read the input file
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Split the file into individual entries
    entries = content.split('----------------------------------------')
    
    for entry in entries:
        if not entry.strip():
            continue
        
        # Extract fields using regular expressions
        gene_match = re.search(r'Gene:\s*(\S+)', entry)
        variation_class_match = re.search(r'Variation Class:\s*(\S+)', entry)
        
        if gene_match and variation_class_match:
            gene = gene_match.group(1)
            variation_class = variation_class_match.group(1)
            gene_data[gene][variation_class] += 1
    
    # Write the summarized data to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['GENE', 'NUM_VARIATIONS', 'VARIATION_TYPES'])
        
        for gene, variations in gene_data.items():
            num_variations = sum(variations.values())
            variation_types = ', '.join([f"{var}({count})" for var, count in variations.items()])
            writer.writerow([gene, num_variations, variation_types])

# Specify input and output file names
input_file = 'data/SNP_data.txt'  # Replace with your input file name
output_file = 'results/SNP_variation_summary.csv'  # Replace with your output file name

# Run the processing function
process_file(input_file, output_file)

print(f"Data has been processed and saved to {output_file}.")

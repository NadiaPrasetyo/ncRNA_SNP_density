import re
import csv
from collections import defaultdict, Counter

# Function to process the text file and extract data
def process_file(input_file, output_file):
    gene_data = defaultdict(lambda: {'variations': Counter(), 'flanks': Counter()})
    
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
        region_type_match = re.search(r'Region Type:\s*(\S+)', entry)
        
        if gene_match and variation_class_match:
            gene = gene_match.group(1)
            variation_class = variation_class_match.group(1)
            region_type = region_type_match.group(1) if region_type_match else 'N/A'
            
            if "flank" in region_type.lower():
                gene_data[gene]['flanks'][variation_class] += 1
            else:
                gene_data[gene]['variations'][variation_class] += 1
    
    # Write the summarized data to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['GENE', 'NUM_VARIATIONS', 'VARIATION_TYPES', 'NUM_FLANKS', 'FLANK_VARIATION_TYPES'])
        
        for gene, data in gene_data.items():
            num_variations = sum(data['variations'].values())
            variation_types = ', '.join([f"{var}({count})" for var, count in data['variations'].items()])
            
            num_flanks = sum(data['flanks'].values())
            flank_variation_types = ', '.join([f"{var}({count})" for var, count in data['flanks'].items()])
            
            writer.writerow([gene, num_variations, variation_types, num_flanks, flank_variation_types])

# Specify input and output file names
# input_file = 'data/SNP_data.txt'  # Replace with your input file name
input_file = 'data/EXTENDED_SNP_data.txt'  
output_file = 'results/EXTENDED_SNP_variation_summary.csv'  # Replace with your output file name

# Run the processing function
process_file(input_file, output_file)

print(f"Data has been processed and saved to {output_file}.")

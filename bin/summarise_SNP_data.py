import re
import csv
from collections import defaultdict, Counter

# Function to process the text file and extract data
def process_file(input_file, output_file):
    gene_data = defaultdict(lambda: {'variations': Counter(), 'flanks': Counter(), 'gene_lengths': [], 'total_lengths': []})
    
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
        original_start_match = re.search(r'Original Start:\s*(\d+)', entry)
        original_end_match = re.search(r'Original End:\s*(\d+)', entry)
        expanded_start_match = re.search(r'Expanded Start:\s*(\d+)', entry)
        expanded_end_match = re.search(r'Expanded End:\s*(\d+)', entry)
        
        # Extract data if matches are found
        if gene_match and variation_class_match and original_start_match and original_end_match and expanded_start_match and expanded_end_match:
            gene = gene_match.group(1)
            variation_class = variation_class_match.group(1)
            region_type = region_type_match.group(1) if region_type_match else 'N/A'
            
            # Extract the start and end positions for both original and expanded ranges
            original_start = int(original_start_match.group(1))
            original_end = int(original_end_match.group(1))
            expanded_start = int(expanded_start_match.group(1))
            expanded_end = int(expanded_end_match.group(1))
            
            # Calculate gene length and total length for the variation
            gene_length = original_end - original_start + 1
            total_length = expanded_end - expanded_start + 1
            
            # Store these lengths for each gene
            gene_data[gene]['gene_lengths'].append(gene_length)
            gene_data[gene]['total_lengths'].append(total_length)
            
            # Count variations and flanks based on region type
            if "flank" in region_type.lower():
                gene_data[gene]['flanks'][variation_class] += 1
            else:
                gene_data[gene]['variations'][variation_class] += 1
    
    # Write the summarized data to the CSV file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['GENE', 'NUM_VARIATIONS', 'VARIATION_TYPES', 'NUM_FLANKS', 'FLANK_VARIATION_TYPES', 'GENE_LENGTH', 'TOTAL_LENGTH'])
        
        for gene, data in gene_data.items():
            num_variations = sum(data['variations'].values())
            variation_types = ', '.join([f"{var}({count})" for var, count in data['variations'].items()])
            
            num_flanks = sum(data['flanks'].values())
            flank_variation_types = ', '.join([f"{var}({count})" for var, count in data['flanks'].items()])
            
            # Calculate the average gene length and total length for the gene
            average_gene_length = sum(data['gene_lengths']) / len(data['gene_lengths']) if data['gene_lengths'] else 0
            average_total_length = sum(data['total_lengths']) / len(data['total_lengths']) if data['total_lengths'] else 0
            
            # Write the summary data with lengths
            writer.writerow([gene, num_variations, variation_types, num_flanks, flank_variation_types, average_gene_length, average_total_length])

# Specify input and output file names
input_file = 'data/EXTENDED_SNP_data.txt'  # Replace with your input file name
output_file = 'results/EXTENDED_SNP_variation_summary.csv'  # Replace with your output file name

# Run the processing function
process_file(input_file, output_file)

print(f"Data has been processed and saved to {output_file}.")

import re
import csv
from collections import defaultdict, Counter

def process_file(input_file, output_file):
    """
    Processes a text file containing gene data, extracts relevant information, and summarizes it in a CSV file.
    
    Args:
        input_file (str): Path to the input text file.
        output_file (str): Path to the output CSV file where results will be saved.
    
    Raises:
        FileNotFoundError: If the input file does not exist.
        IOError: If there is an error reading from the input file or writing to the output file.
        ValueError: If the input file contains invalid data or cannot be parsed correctly.
    """
    # Dictionary to store gene data with default structures for variations, flanks, and lengths
    gene_data = defaultdict(lambda: {
        'variations': Counter(),
        'flanks': Counter(),
        'gene_lengths': [],
        'total_lengths': []
    })
    
    try:
        # Read the input file
        with open(input_file, 'r') as f:
            content = f.read()
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{input_file}' does not exist.")
    except IOError as e:
        raise IOError(f"Error reading from '{input_file}': {e}")
    
    # Split the file into individual entries
    entries = content.split('----------------------------------------')
    
    for entry in entries:
        if not entry.strip():
            continue  # Skip empty entries
        
        try:
            # Extract fields using regular expressions
            gene_match = re.search(r'Gene:\s*(\S+)', entry)
            variation_class_match = re.search(r'Variation Class:\s*(\S+)', entry)
            region_type_match = re.search(r'Region Type:\s*(.*)', entry)  # Match full region type
            original_start_match = re.search(r'Original Start:\s*(\d+)', entry)
            original_end_match = re.search(r'Original End:\s*(\d+)', entry)
            expanded_start_match = re.search(r'Expanded Start:\s*(\d+)', entry)
            expanded_end_match = re.search(r'Expanded End:\s*(\d+)', entry)
            
            # Check if all necessary fields are present
            if not (gene_match and variation_class_match and original_start_match and original_end_match and expanded_start_match and expanded_end_match):
                raise ValueError("One or more required fields are missing in an entry.")
            
            # Extract data
            gene = gene_match.group(1)
            variation_class = variation_class_match.group(1)
            region_type = region_type_match.group(1).strip() if region_type_match else 'N/A'
            region_type = region_type.replace(" ", "")  # Removes spaces between words
            
            # Parse and calculate lengths
            original_start = int(original_start_match.group(1))
            original_end = int(original_end_match.group(1))
            expanded_start = int(expanded_start_match.group(1))
            expanded_end = int(expanded_end_match.group(1))
            gene_length = original_end - original_start + 1
            total_length = expanded_end - expanded_start + 1
            
            # Update gene data
            gene_data[gene]['gene_lengths'].append(gene_length)
            gene_data[gene]['total_lengths'].append(total_length)
            
            # Count variations and flanks
            if "flank" in region_type.lower():
                gene_data[gene]['flanks'][variation_class] += 1
            else:
                gene_data[gene]['variations'][variation_class] += 1
        
        except ValueError as e:
            print(f"Error processing entry: {e}")
            continue  # Skip the problematic entry
    
    try:
        # Write summarized data to the CSV file
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['GENE', 'NUM_VARIATIONS', 'VARIATION_TYPES', 'NUM_FLANKS', 'FLANK_VARIATION_TYPES', 'GENE_LENGTH', 'TOTAL_LENGTH'])
            
            for gene, data in gene_data.items():
                num_variations = sum(data['variations'].values())
                variation_types = ', '.join([f"{var}({count})" for var, count in data['variations'].items()])
                num_flanks = sum(data['flanks'].values())
                flank_variation_types = ', '.join([f"{var}({count})" for var, count in data['flanks'].items()])
                
                # Calculate average lengths
                average_gene_length = (sum(data['gene_lengths']) / len(data['gene_lengths'])) if data['gene_lengths'] else 0
                average_total_length = (sum(data['total_lengths']) / len(data['total_lengths'])) if data['total_lengths'] else 0
                
                writer.writerow([gene, num_variations, variation_types, num_flanks, flank_variation_types, average_gene_length, average_total_length])
    except IOError as e:
        raise IOError(f"Error writing to '{output_file}': {e}")
    
    print(f"Data has been processed and saved to {output_file}.")

# Specify input and output file names
input_file = 'data/EXTENDED_SNP_data.txt'  # Replace with your input file name
output_file = 'results/EXTENDED_SNP_variation_summary.csv'  # Replace with your output file name

# Run the processing function
try:
    process_file(input_file, output_file)
except Exception as e:
    print(f"An error occurred: {e}")

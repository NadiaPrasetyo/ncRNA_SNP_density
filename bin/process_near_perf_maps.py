import re
import csv
from collections import defaultdict

def extract_chromosome_from_filename(filename):
    # Use a regular expression to extract '-chr#' from the filename
    match = re.search(r'-chr(\d+)', filename)
    if match:
        return f"chr{match.group(1)}"  # Return the chromosome in 'chrX' format
    return None  # Return None if no match found

def process_sam_file(sam_file):
    # Store results: gene name -> {read_length, unique chromosomes for the hits, non-perfect hits, perfect hit count}
    results = defaultdict(lambda: {'read_length': None, 'chromosomes': set(), 'non_perfect_hits': 0, 'perfect_hit_count': 0})
    
    # Read the .sam file
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Skip header lines
                continue

            # Split the line into columns (assuming column 0 is the filename, and the rest are standard SAM fields)
            columns = line.strip().split('\t')
            
            # Ensure the line has at least 12 columns (standard SAM format)
            if len(columns) < 12:
                continue  # Skip lines that don't have the expected number of columns

            # Extract necessary fields from the modified columns
            filename = columns[0]  # Filename is in the first column
            gene_name = columns[1]  # Gene name (column 1)
            chromosome = extract_chromosome_from_filename(filename)
            position = columns[3]  # Position (column 3)
            cigar = columns[6]  # CIGAR string (column 6)
            sequence = columns[9]  # Sequence (column 9)
            quality_score = int(columns[5])  # Quality score (column 5)

            # Initialize values
            alignment_score = 0
            num_mismatches = 0

            # Try to extract the AS:i:<value> and NM:i:<value> fields (they start from column 11 and onward)
            for field in columns[11:]:
                if field.startswith('AS:i:'):
                    try:
                        alignment_score = int(field.split(':')[2])  # AS:i:<value>
                    except ValueError:
                        continue  # Skip malformed AS:i: fields
                elif field.startswith('NM:i:'):
                    try:
                        num_mismatches = int(field.split(':')[2])  # NM:i:<value>
                    except ValueError:
                        continue  # Skip malformed NM:i: fields

            # Check if the hit is a perfect hit (no 'X' for mismatches, quality score 42)
            if 'X' not in cigar and num_mismatches == 0 and quality_score == 42:
                # Get the read length from the CIGAR string
                read_length = sum(int(x) for x in re.findall(r'(\d+)', cigar))  # Length from CIGAR
                
                # For the first perfect hit, store the read length
                if results[gene_name]['read_length'] is None:
                    results[gene_name]['read_length'] = read_length
                
                # Store the unique chromosome where the perfect hit occurred
                if chromosome not in results[gene_name]['chromosomes']:
                    results[gene_name]['chromosomes'].add(chromosome)
                    # Increment the perfect hit count for this gene
                    results[gene_name]['perfect_hit_count'] += 1
            else:
                # Non-perfect hit
                results[gene_name]['non_perfect_hits'] += 1
    
    return results

def generate_summary(results):
    # Generate a summary output
    summary = []
    
    for gene_name, data in results.items():
        # Get read length, unique chromosomes, and non-perfect hit count
        read_length = data['read_length']
        chromosomes = sorted(data['chromosomes'])
        non_perfect_hits = data['non_perfect_hits']
        perfect_hit_count = data['perfect_hit_count']
        
        # Format the summary for each gene
        if read_length:
            # Create a string for the chromosomes without repeating the read length
            chromosomes_summary = ", ".join([f"{chr}" for chr in chromosomes])
            summary.append([gene_name, f"{perfect_hit_count} ({read_length}bp, {chromosomes_summary})", non_perfect_hits])
    
    return summary

def write_summary_to_csv(summary, output_file):
    # Write the summary to a CSV file
    with open(output_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        # CSV header
        writer.writerow(['Gene Name', 'Number of Unique Mapping', 'No. of Near Perfect Matches'])
        # Write each row of the summary
        writer.writerows(summary)

# Example Usage
sam_file = 'data/aligned_reads.sam'  # Provide your .sam file path
results = process_sam_file(sam_file)
summary = generate_summary(results)

# Write the summary of results to a CSV file
output_file = 'results/GED-MAP-alignment-summarised.csv'  # Specify the desired output file name
write_summary_to_csv(summary, output_file)

print(f"Summary has been written to {output_file}")

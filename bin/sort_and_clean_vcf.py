import os

def is_valid_chromosome(chromosome):
    """
    Check if the chromosome is a valid numeric chromosome (1-22).
    
    Args:
        chromosome (str): Chromosome string (e.g., "chr1", "chrX").
    
    Returns:
        bool: True if chromosome is between "chr1" and "chr22", False otherwise.
    """
    chrom_name = chromosome[3:]  # Skip the "chr" prefix
    if chrom_name.isdigit():
        return 1 <= int(chrom_name) <= 22  # Only numeric chromosomes 1-22 are valid
    return False

def parse_chromosome(chromosome):
    """
    Parse chromosome for sorting purposes.
    
    Args:
        chromosome (str): Chromosome string (e.g., "chr1", "chrX").
    
    Returns:
        int or str: Numeric part of the chromosome for valid chromosomes, else the original string.
    """
    chrom_name = chromosome[3:]
    return int(chrom_name) if chrom_name.isdigit() else chromosome

def process_vcf(input_file, output_file):
    """
    Process a VCF file by filtering and sorting its entries based on specific criteria.
    
    Args:
        input_file (str): Path to the input VCF file.
        output_file (str): Path to save the filtered and sorted VCF file.
    """
    try:
        # Read the input VCF file
        with open(input_file, 'r') as infile:
            lines = infile.readlines()
        
        # Separate header and data lines
        header_lines = [line for line in lines if line.startswith('#')]
        data_lines = [line for line in lines if not line.startswith('#')]
        
        # Filter data lines based on REF and ALT lengths and chromosome validity
        filtered_data_lines = [
            line for line in data_lines
            if len(line.split('\t')[3]) <= 40  # REF length
            and len(line.split('\t')[4]) <= 40  # ALT length
            and is_valid_chromosome(line.split('\t')[0])  # Chromosome validity
        ]
        
        # Sort filtered data lines by chromosome number
        sorted_data_lines = sorted(filtered_data_lines, key=lambda x: parse_chromosome(x.split('\t')[0]))
        
        # Write the results to the output file, including the original header
        with open(output_file, 'w') as outfile:
            outfile.writelines(header_lines)
            outfile.writelines(sorted_data_lines)
        
        print(f"Sorted and filtered VCF file has been saved to {output_file}")
    
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except IOError as e:
        print(f"Error reading or writing file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# File paths
input_file = 'data/decomposed.vcf'
output_file = 'data/sorted_all_variants.vcf'

# Ensure the input file exists before processing
if os.path.exists(input_file):
    process_vcf(input_file, output_file)
else:
    print(f"Input file '{input_file}' does not exist. Please provide a valid file path.")

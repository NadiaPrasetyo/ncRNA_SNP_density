import os

def is_valid_chromosome(chromosome, validate=True):
    """
    Check if the chromosome is a valid numeric chromosome (1-22).
    
    Args:
        chromosome (str): Chromosome string (e.g., "chr1", "chrX").
        validate (bool): Whether to validate chromosomes.
    
    Returns:
        bool: True if chromosome is valid or validation is disabled, False otherwise.
    """
    if not validate:
        return True  # Skip validation if disabled
    
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
        tuple: A tuple for sorting. Numeric chromosomes are sorted numerically,
               and others are sorted lexicographically after the numeric ones.
    """
    chrom_name = chromosome[3:]  # Skip the "chr" prefix
    if chrom_name.isdigit():
        return (0, int(chrom_name))  # Numeric chromosomes sort first
    return (1, chrom_name)  # Non-numeric chromosomes sort afterward


def process_vcf(input_file, output_file, validate_chromosome=True):
    """
    Process a VCF file by filtering and sorting its entries based on specific criteria.
    
    Args:
        input_file (str): Path to the input VCF file.
        output_file (str): Path to save the filtered and sorted VCF file.
        validate_chromosome (bool): Whether to validate chromosomes for filtering.
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
            and is_valid_chromosome(line.split('\t')[0], validate=validate_chromosome)  # Chromosome validity
        ]
        
        # Sort filtered data lines by chromosome number
        sorted_data_lines = sorted(filtered_data_lines, key=lambda x: parse_chromosome(x.split('\t')[0]))
        
        # Write the results to the output file, including the original header
        with open(output_file, 'w') as outfile:
            outfile.writelines(header_lines)
            outfile.writelines(sorted_data_lines)
        
        print(f"Processed '{input_file}' -> '{output_file}'")
    
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
    except IOError as e:
        print(f"Error reading or writing file: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def process_vcf_folder(input_folder, output_folder, validate_chromosome=True):
    """
    Process all VCF files in a given folder.
    
    Args:
        input_folder (str): Path to the folder containing VCF files.
        output_folder (str): Path to the folder where processed files will be saved.
        validate_chromosome (bool): Whether to validate chromosomes for filtering.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # Create output folder if it doesn't exist
    
    # Iterate through files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith('.vcf'):  # Process only VCF files
            input_file = os.path.join(input_folder, file_name)
            output_file = os.path.join(output_folder, f"processed_{file_name}")
            process_vcf(input_file, output_file, validate_chromosome=validate_chromosome)

def main():
    # User configuration
    single_file_mode = False  # Set to True to process a single file, False to process a folder
    validate_chromosome = False  # Set to False to skip chromosome validation

    if single_file_mode:
        # Single file processing
        input_file = 'data/decomposed.vcf'
        output_file = 'data/sorted_all_variants.vcf'
        
        if os.path.exists(input_file):
            process_vcf(input_file, output_file, validate_chromosome=validate_chromosome)
        else:
            print(f"Input file '{input_file}' does not exist. Please provide a valid file path.")
    else:
        # Folder processing
        input_folder = 'data/VCF'
        output_folder = 'data/VCF/processed'
        
        if os.path.exists(input_folder):
            process_vcf_folder(input_folder, output_folder, validate_chromosome=validate_chromosome)
        else:
            print(f"Input folder '{input_folder}' does not exist. Please provide a valid folder path.")

if __name__ == "__main__":
    main()

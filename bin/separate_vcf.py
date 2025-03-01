import os

def split_vcf_by_chrom(input_vcf=None, input_folder=None):
    """
    Splits VCF (Variant Call Format) files into separate files for each chromosome.
    
    Parameters:
    input_vcf (str): Path to a single VCF file (optional, if input_folder is provided, this will be ignored).
    input_folder (str): Path to the folder containing VCF files (optional, if input_vcf is provided, this will be ignored).

    Output:
    Creates a separate VCF file for each chromosome in an output directory named 'data/VCF/'.
    """
    if not input_vcf and not input_folder:
        print("Error: You must provide either a file or a folder as input.")
        return

    # Ensure the output directory exists
    output_dir = "data/VCF/"
    os.makedirs(output_dir, exist_ok=True)

    # Handle the case where a single file is provided
    if input_vcf:
        if os.path.isfile(input_vcf):
            print(f"Processing single file: {input_vcf}")
            process_vcf(input_vcf, output_dir)
        else:
            print(f"Error: The file '{input_vcf}' was not found.")
    
    # Handle the case where a folder is provided
    if input_folder:
        if os.path.isdir(input_folder):
            print(f"Processing files in folder: {input_folder}")
            for file_name in os.listdir(input_folder):
                if file_name.endswith(".vcf"):
                    input_vcf_path = os.path.join(input_folder, file_name)
                    print(f"Processing file: {input_vcf_path}")
                    process_vcf(input_vcf_path, output_dir)
        else:
            print(f"Error: The folder '{input_folder}' was not found.")


def process_vcf(input_vcf_path, output_dir):
    """
    Processes a VCF file and splits it by chromosome.

    Parameters:
    input_vcf_path (str): Path to the VCF file to be processed.
    output_dir (str): Path to the output directory where split VCF files will be saved.
    """
    chrom_files = {}

    try:
        with open(input_vcf_path, 'r') as infile:
            header_lines = []
            for line in infile:
                if line.startswith("##"):
                    # Collect meta-information header lines
                    header_lines.append(line)
                elif line.startswith("#CHROM"):
                    # This is the column header; collect it to write to all chromosome files
                    header_lines.append(line)
                else:
                    # Process each data line
                    try:
                        fields = line.strip().split("\t")
                        chrom = fields[0]  # The chromosome column is the first column (index 0)
                        
                        # If chrom is not in the dictionary, create a new file for it
                        if chrom not in chrom_files:
                            chrom_file_path = os.path.join(output_dir, f"{os.path.basename(input_vcf_path)}_{chrom}.vcf")
                            chrom_files[chrom] = open(chrom_file_path, 'w')
                            
                            # Write the header lines to the new file
                            for header in header_lines:
                                chrom_files[chrom].write(header)
                        
                        # Write the variant line to the corresponding file
                        chrom_files[chrom].write(line)
                    
                    except IndexError:
                        print(f"Skipping malformed line: {line.strip()}")
                    except Exception as e:
                        print(f"Error processing line: {line.strip()}. Error: {e}")
        
        print(f"Finished processing file: {input_vcf_path}")

    except FileNotFoundError:
        print(f"Error: The input file '{input_vcf_path}' was not found.")
    except PermissionError:
        print(f"Error: Permission denied while accessing '{input_vcf_path}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    finally:
        # Ensure all file handlers are closed
        for chrom, file in chrom_files.items():
            try:
                file.close()
            except Exception as e:
                print(f"Error closing file for chromosome '{chrom}': {e}")


# Example usage:
# For a single VCF file input:
# split_vcf_by_chrom(input_vcf="data/sorted_all_variants.vcf")

# For a folder containing VCF files:
split_vcf_by_chrom(input_folder="data/VCF")

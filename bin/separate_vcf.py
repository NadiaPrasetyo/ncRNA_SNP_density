import os

def split_vcf_by_chrom(input_vcf):
    """
    Splits a VCF (Variant Call Format) file into separate files for each chromosome or gene.

    Parameters:
    input_vcf (str): Path to the input VCF file.

    Output:
    Creates a separate VCF file for each chromosome in the specified output directory.
    """
    # Dictionary to hold file handlers for each chromosome
    chrom_files = {}
    
    try:
        # Ensure the output directory exists
        output_dir = "data/VCF/"
        os.makedirs(output_dir, exist_ok=True)

        with open(input_vcf, 'r') as infile:
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
                            chrom_file_path = os.path.join(output_dir, f"{chrom}.vcf")
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
        
        print("VCF splitting completed successfully.")

    except FileNotFoundError:
        print(f"Error: The input file '{input_vcf}' was not found.")
    except PermissionError:
        print(f"Error: Permission denied while accessing '{input_vcf}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    
    finally:
        # Ensure all file handlers are closed
        for chrom, file in chrom_files.items():
            try:
                file.close()
            except Exception as e:
                print(f"Error closing file for chromosome '{chrom}': {e}")

# Input VCF file
input_vcf_path = "data/sorted_all_variants.vcf"  # Replace with your VCF file path

# Call the function to split the VCF by chromosome
split_vcf_by_chrom(input_vcf_path)

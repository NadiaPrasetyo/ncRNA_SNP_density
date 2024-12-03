import os

def split_vcf_by_chrom(input_vcf):
    # Create a dictionary to hold file handlers for each gene/chromosome
    chrom_files = {}

    with open(input_vcf, 'r') as infile:
        header_lines = []
        for line in infile:
            if line.startswith("##"):
                # Collect header lines
                header_lines.append(line)
            elif line.startswith("#CHROM"):
                # This is the column header, we need to process it and write to each file
                header_lines.append(line)
            else:
                # Process each data line
                fields = line.strip().split("\t")
                chrom = fields[0]  # CHROM column is the 1st column (index 0)
                
                # If chrom is not in the dictionary, create a new file for it
                if chrom not in chrom_files:
                    chrom_files[chrom] = open(f"data/VCF/{chrom}_variants.vcf", 'w')
                    # Write the header lines to the new file
                    for header in header_lines:
                        chrom_files[chrom].write(header)
                
                # Write the variant line to the corresponding file
                chrom_files[chrom].write(line)

    # Close all file handlers
    for chrom_file in chrom_files.values():
        chrom_file.close()


# Input VCF file
input_vcf_path = "data/CLEAN_filtered_variants.vcf"  # Replace with your VCF file path

# Split the VCF by chromosome (or gene)
split_vcf_by_chrom(input_vcf_path)

print("VCF splitting complete. Check the output files.")

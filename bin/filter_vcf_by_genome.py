import pandas as pd

def filter_vcf_by_genome(input_vcf, output_dir):
    """
    Filters a VCF file to create individual VCF files for each genome,
    trimming columns to include only the necessary fields and adjusting ALT alleles.

    Args:
        input_vcf (str): Path to the input VCF file.
        output_dir (str): Directory where the filtered VCF files will be saved.
    """
    # Read the VCF file, separating header and data
    header_lines = []
    with open(input_vcf, 'r') as file:
        for line in file:
            if line.startswith("##"):
                header_lines.append(line)
            elif line.startswith("#"):
                header_line = line.strip()  # This contains the column headers
                break

    # Load the data using the extracted header
    column_names = header_line.lstrip("#").split("\t")
    vcf_data = pd.read_csv(input_vcf, sep='\t', comment='#', names=column_names, skiprows=len(header_lines))

    # Debugging: Verify column names
    print("Parsed Columns:", vcf_data.columns)

    # Ensure column names are stripped of extra whitespace
    vcf_data.columns = vcf_data.columns.str.strip()

    # Extract the list of genomes
    genomes = vcf_data.columns[10:]  # All columns after FORMAT

    # Process each genome
    for genome in genomes:
        # Ensure all values are treated as strings
        genome_data = vcf_data[genome].astype(str)
        
        # Check if the genome has the variation
        genome_variations = vcf_data[genome_data.apply(
            lambda x: any(
                allele not in {'0', '.', '|'} for allele in x.split('|')
            )
        )]
        
        # Trim to necessary columns
        trimmed_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'CHM13', genome]
        genome_variations = genome_variations[trimmed_columns]

        # Adjust ALT field based on the genome's alleles
        def adjust_alt(row):
            # Split ALT into a list of alleles
            alt_alleles = row['ALT'].split(',')
            # Get the genotype field for the genome
            genotype = row[genome]
            # Extract relevant allele indices from the genotype
            relevant_indices = [
                int(allele) - 1
                for allele in genotype.replace('|', '/').split('/')
                if allele.isdigit() and int(allele) > 0
            ]
            # Filter ALT alleles by relevant indices
            filtered_alt = [alt_alleles[i] for i in relevant_indices]
            return ','.join(filtered_alt) if filtered_alt else '.'

        # Apply ALT adjustment
        genome_variations['ALT'] = genome_variations.apply(adjust_alt, axis=1)

        # Prepare output VCF
        genome_header = header_lines[:]
        genome_header.append(f"##Filtered for genome {genome}\n")
        
        # Save the filtered data
        output_file = f"{output_dir}/{genome}.vcf"
        with open(output_file, 'w') as out_vcf:
            out_vcf.writelines(genome_header)
            genome_variations.to_csv(out_vcf, sep='\t', index=False, header=True)

# Usage
filter_vcf_by_genome("data/filtered_variants.vcf", "data/VCF")

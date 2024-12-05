import csv
import re

# Increase the field size limit to a high value
csv.field_size_limit(10**7)

# Input and output file paths
input_file = 'data/EXTENDED_filtered_variants.csv'  # Replace with your actual input file path
output_file = 'data/EXTENDED_pangenome_summary.csv'

# Function to extract TYPE from INFO field
def extract_type(info):
    match = re.search(r"TYPE=([^;]+)", info)
    return match.group(1) if match else ""

# Process the CSV file
with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile)
    
    # Extract all genome columns
    genome_columns = [col for col in reader.fieldnames if col not in 
                      ('GENE', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                       'ORIGINAL_START', 'ORIGINAL_END', 'EXTENDED_START', 'EXTENDED_END')]  # Exclude new fields
    
    # Prepare the writer for the output file, adding new fields for start and end positions
    fieldnames = ['GENE', 'CHROM', 'POS', 'TYPE', 'SNP_TYPE', 'GENOME_DETAILS', 
                  'ORIGINAL_START', 'ORIGINAL_END', 'EXTENDED_START', 'EXTENDED_END']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()

    # Process each row
    for row in reader:
        gene = row['GENE']
        chrom = row['CHROM']
        pos = row['POS']
        ref = row['REF']
        alt = row['ALT'].split(",")  # Handle multiple alternative alleles
        info = row['INFO']
        
        # Extract TYPE from the INFO field
        snp_type_info = extract_type(info)
        
        # If TYPE is not found, assign "snp" as default
        if not snp_type_info:
            snp_type_info = "snp"

        # Construct SNP_TYPE (e.g., A>G, A>T)
        snp_types = [f"{ref}>{a}" for a in alt]

        # Gather genome details
        genome_details = []
        for genome in genome_columns:
            gt = row[genome]
            if gt == '0|0' or gt == '0':  # No SNP
                continue
            elif gt in ('1|1', '2|2', '3|3'):  # Homozygous
                allele_index = int(gt.split("|")[0]) - 1  # 1-based index for ALT
                if 0 <= allele_index < len(snp_types):
                    genome_details.append(f"{genome} (Homozygous alt {allele_index + 1})")
            elif re.match(r"0\|[1-9]|[1-9]\|0", gt):  # Heterozygous
                allele_index = max(int(gt.split("|")[0]), int(gt.split("|")[1])) - 1
                if 0 <= allele_index < len(snp_types):
                    genome_details.append(f"{genome} (Heterozygous alt {allele_index + 1})")
        
        # Write to the output file, including the start/end positions
        writer.writerow({
            'GENE': gene,
            'CHROM': chrom,
            'POS': pos,
            'TYPE': snp_type_info,
            'SNP_TYPE': ",".join(snp_types),
            'GENOME_DETAILS': "; ".join(genome_details),
            'ORIGINAL_START': row['ORIGINAL_START'],
            'ORIGINAL_END': row['ORIGINAL_END'],
            'EXTENDED_START': row['EXTENDED_START'],
            'EXTENDED_END': row['EXTENDED_END']
        })

print(f"Processing complete. The output file is saved as {output_file}.")

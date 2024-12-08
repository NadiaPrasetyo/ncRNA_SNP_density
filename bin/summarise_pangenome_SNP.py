import csv
import re

# Increase the field size limit to handle large CSV files
csv.field_size_limit(10**7)

# Input and output file paths
input_file = 'data/EXTENDED_filtered_variants.csv'  # Replace with your actual input file path
output_file = 'data/EXTENDED_pangenome_summary.csv'

def extract_type(info):
    """
    Extracts the TYPE value from the INFO field of a VCF-like entry.
    
    Args:
        info (str): The INFO field containing metadata about the variant.
    
    Returns:
        str: The extracted TYPE value, or an empty string if not found.
    """
    match = re.search(r"TYPE=([^;]+)", info)
    return match.group(1) if match else ""

try:
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        try:
            reader = csv.DictReader(infile)
        except csv.Error as e:
            raise ValueError(f"Error reading CSV input file: {e}")

        # Extract all genome columns, excluding known metadata columns
        genome_columns = [col for col in reader.fieldnames if col not in 
                          ('GENE', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                           'ORIGINAL_START', 'ORIGINAL_END', 'EXTENDED_START', 'EXTENDED_END')]
        
        # Prepare the writer for the output file
        fieldnames = ['GENE', 'CHROM', 'POS', 'TYPE', 'SNP_TYPE', 'GENOME_DETAILS', 
                      'ORIGINAL_START', 'ORIGINAL_END', 'EXTENDED_START', 'EXTENDED_END']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            try:
                # Extract relevant fields
                gene = row['GENE']
                chrom = row['CHROM']
                pos = row['POS']
                ref = row['REF']
                alt = row['ALT'].split(",")  # Handle multiple alternative alleles
                info = row['INFO']
                
                # Extract TYPE from the INFO field
                snp_type_info = extract_type(info)
                if not snp_type_info:
                    snp_type_info = "snp"

                # Construct SNP_TYPE (e.g., A>G, A>T)
                snp_types = [f"{ref}>{a}" for a in alt]

                # Gather genome details
                genome_details = []
                for genome in genome_columns:
                    gt = row.get(genome, "")
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
                
                # Write processed data to the output file
                writer.writerow({
                    'GENE': gene,
                    'CHROM': chrom,
                    'POS': pos,
                    'TYPE': snp_type_info,
                    'SNP_TYPE': ",".join(snp_types),
                    'GENOME_DETAILS': "; ".join(genome_details),
                    'ORIGINAL_START': row.get('ORIGINAL_START', ''),
                    'ORIGINAL_END': row.get('ORIGINAL_END', ''),
                    'EXTENDED_START': row.get('EXTENDED_START', ''),
                    'EXTENDED_END': row.get('EXTENDED_END', '')
                })
            except KeyError as e:
                print(f"Skipping row due to missing key: {e}")
            except Exception as e:
                print(f"Error processing row {row}: {e}")
    
    print(f"Processing complete. The output file is saved as {output_file}.")
except FileNotFoundError:
    print(f"The input file '{input_file}' was not found.")
except IOError as e:
    print(f"Error accessing files: {e}")
except Exception as e:
    print(f"An unexpected error occurred: {e}")

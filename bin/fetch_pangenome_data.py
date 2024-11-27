import csv

# Define the file path the csv file of genes and locations with high SNP density
csv_file_path = r"data\SNP-densities-and-RNA.csv"

# Define the file path to your VCF file
vcf_file_path = r"data\decomposed.vcf"

# Define the output file path for the filtered variants
output_csv_file_path = "results/filtered_variants.csv"

# Column headers for the output CSV
output_headers = [
    "GENE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    "CHM13", "HG00438", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741",
    "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258",
    "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055",
    "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559",
    "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886",
    "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579",
    "NA18906", "NA20129", "NA21309"
]

# Read the chromosome, start, and end ranges from the CSV file
def read_gene_ranges(csv_file_path):
    gene_ranges = []
    try:
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                gene_ranges.append({
                    "chromosome": row["Chromosome"],
                    "start": int(row["Start"]),
                    "end": int(row["End"]),
                    "gene_name": row["GeneName"]
                })
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file_path}")
        return []
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return []
    return gene_ranges

# Filter VCF lines based on gene ranges and prepare them for CSV output
def filter_vcf_and_prepare_csv(vcf_file_path, gene_ranges, lines_to_skip=2595):
    csv_rows = []
    try:
        with open(vcf_file_path, 'r') as vcf_file:
            # Skip the first 2595 lines
            for _ in range(lines_to_skip):
                vcf_file.readline()
            
            # Process the remaining lines
            for line in vcf_file:
                if line.startswith("#"):  # Skip header lines
                    continue
                fields = line.strip().split("\t")
                
                # Extract VCF fields
                chrom = fields[0]
                pos = int(fields[1])
                vcf_data = {
                    "CHROM": chrom,
                    "POS": pos,
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "QUAL": fields[5],
                    "FILTER": fields[6],
                    "INFO": fields[7],
                    "FORMAT": fields[8],
                    "GENOTYPES": fields[9:]  # Remaining columns are genotype data
                }

                # Check if this variant falls within any gene range
                for gene in gene_ranges:
                    if chrom == gene["chromosome"] and gene["start"] <= pos <= gene["end"]:
                        # Print progress
                        print(f"Processing gene: {gene['gene_name']}, chromosome location: {gene['chromosome']}:{gene['start']}-{gene['end']}")
                        print(f"Variation found: ID={vcf_data['ID']}, POS={vcf_data['POS']}")
                        
                        # Add a new row to the CSV
                        csv_row = {
                            "GENE": gene["gene_name"],
                            "CHROM": vcf_data["CHROM"],
                            "POS": vcf_data["POS"],
                            "ID": vcf_data["ID"],
                            "REF": vcf_data["REF"],
                            "ALT": vcf_data["ALT"],
                            "QUAL": vcf_data["QUAL"],
                            "FILTER": vcf_data["FILTER"],
                            "INFO": vcf_data["INFO"],
                            "FORMAT": vcf_data["FORMAT"]
                        }
                        # Add genotypes to the row
                        for idx, sample in enumerate(output_headers[10:]):  # Genotype headers start at index 10
                            csv_row[sample] = vcf_data["GENOTYPES"][idx] if idx < len(vcf_data["GENOTYPES"]) else "."
                        csv_rows.append(csv_row)
                        break  # Stop checking other ranges once a match is found
    except FileNotFoundError:
        print(f"VCF file not found at {vcf_file_path}")
    except Exception as e:
        print(f"Error reading VCF file: {e}")
    return csv_rows

# Main workflow
def main():
    # Read gene ranges from the CSV file
    gene_ranges = read_gene_ranges(csv_file_path)
    if not gene_ranges:
        print("No gene ranges available to process.")
        return

    # Filter the VCF file and prepare rows for the CSV
    csv_rows = filter_vcf_and_prepare_csv(vcf_file_path, gene_ranges)

    # Write the filtered variants to a new CSV file
    try:
        with open(output_csv_file_path, 'w', newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=output_headers)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Filtered variants written to {output_csv_file_path}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

if __name__ == "__main__":
    main()
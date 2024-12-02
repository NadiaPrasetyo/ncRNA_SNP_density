import csv

# Define file paths
csv_file_path = r"data/SNP-densities-and-RNA.csv"
vcf_file_path = r"data/decomposed.vcf"
output_csv_file_path = "results/filtered_variants.csv"
output_vcf_file_path = "data/filtered_variants.vcf"

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

# Read gene ranges from the CSV file
def read_and_sort_gene_ranges(csv_file_path):
    print("Reading and sorting gene ranges from the CSV file...")
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
        # Sort by chromosome and start position
        gene_ranges.sort(key=lambda x: (x["chromosome"], x["start"]))
        print(f"Successfully sorted {len(gene_ranges)} gene ranges.")
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file_path}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return gene_ranges

# Filter VCF file based on sorted ranges
def filter_vcf_by_ranges(vcf_file_path, gene_ranges, lines_to_skip=2595):
    print("Filtering VCF file based on ranges...")
    csv_rows = []
    vcf_lines = []
    header_lines = []
    current_range_idx = 0
    total_variants_matched = 0

    try:
        with open(vcf_file_path, 'r') as vcf_file:
            # Skip the first 2595 lines
            for i in range(lines_to_skip):
                line = vcf_file.readline()
                if line.startswith("#"):  # Retain header lines for VCF output
                    header_lines.append(line)

            # Process VCF file
            for line in vcf_file:
                if line.startswith("#"):
                    continue

                # Parse VCF line
                fields = line.strip().split("\t")
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
                    "GENOTYPES": fields[9:]
                }

                # Match against gene ranges
                while current_range_idx < len(gene_ranges):
                    gene = gene_ranges[current_range_idx]

                    # If the variant is before the current range, skip it
                    if chrom < gene["chromosome"] or (chrom == gene["chromosome"] and pos < gene["start"]):
                        break

                    # If the variant is beyond the current range, move to the next range
                    if chrom > gene["chromosome"] or (chrom == gene["chromosome"] and pos > gene["end"]):
                        current_range_idx += 1
                        continue

                    # If the variant is within the range, save it
                    if chrom == gene["chromosome"] and gene["start"] <= pos <= gene["end"]:
                        print(f"Match found: Gene={gene['gene_name']}, Location={chrom}:{pos}")
                        total_variants_matched += 1

                        # Add data for CSV output
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
                        for idx, sample in enumerate(output_headers[10:]):
                            csv_row[sample] = vcf_data["GENOTYPES"][idx] if idx < len(vcf_data["GENOTYPES"]) else "."
                        csv_rows.append(csv_row)

                        # Add data for VCF output
                        vcf_lines.append(line.strip())
                        break
    except FileNotFoundError:
        print(f"VCF file not found at {vcf_file_path}")
    except Exception as e:
        print(f"Error reading VCF file: {e}")

    print(f"Total variants matched: {total_variants_matched}")
    return header_lines, csv_rows, vcf_lines

# Write filtered variants to CSV
def write_csv(csv_rows):
    print("Writing the filtered variants to the output CSV file...")
    try:
        with open(output_csv_file_path, 'w', newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=output_headers)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Filtered variants successfully written to {output_csv_file_path}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

# Write filtered variants to VCF
def write_vcf(header_lines, vcf_lines):
    print("Writing the filtered variants to the output VCF file...")
    try:
        with open(output_vcf_file_path, 'w') as output_file:
            output_file.writelines(header_lines)
            output_file.write("\n".join(vcf_lines) + "\n")
        print(f"Filtered variants successfully written to {output_vcf_file_path}")
    except Exception as e:
        print(f"Error writing VCF file: {e}")

# Main workflow
def main():
    print("Starting the filtering workflow...")

    # Step 1: Read and sort gene ranges
    gene_ranges = read_and_sort_gene_ranges(csv_file_path)
    if not gene_ranges:
        print("No gene ranges available to process. Exiting.")
        return

    # Step 2: Filter the VCF file
    header_lines, csv_rows, vcf_lines = filter_vcf_by_ranges(vcf_file_path, gene_ranges)

    # Step 3: Choose output format
    output_format = input("Enter output format (csv/vcf): ").strip().lower()
    if output_format == "csv":
        write_csv(csv_rows)
    elif output_format == "vcf":
        write_vcf(header_lines, vcf_lines)
    else:
        print("Invalid output format. Exiting.")

    print("Workflow completed.")

if __name__ == "__main__":
    main()

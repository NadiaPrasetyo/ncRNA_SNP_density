import csv
import re

# Define file paths
csv_file_path = r"data/SNP-densities-and-RNA.csv"
vcf_file_path = r"data/decomposed.vcf"
output_csv_file_path = "data/EXTENDED_filtered_variants.csv"
output_vcf_file_path = "data/EXTENDED_filtered_variants.vcf"

# Column headers for the output CSV (add new columns for start and end positions)
output_headers = [
    "GENE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",  "ORIGINAL_START", "ORIGINAL_END", 
    "EXTENDED_START", "EXTENDED_END", "HG00438", "HG00621", "HG00673", "HG00733", "HG00735", "HG00741",
    "HG01071", "HG01106", "HG01109", "HG01123", "HG01175", "HG01243", "HG01258",
    "HG01358", "HG01361", "HG01891", "HG01928", "HG01952", "HG01978", "HG02055",
    "HG02080", "HG02109", "HG02145", "HG02148", "HG02257", "HG02486", "HG02559",
    "HG02572", "HG02622", "HG02630", "HG02717", "HG02723", "HG02818", "HG02886",
    "HG03098", "HG03453", "HG03486", "HG03492", "HG03516", "HG03540", "HG03579",
    "NA18906", "NA20129", "NA21309"
]

def sort_chromosomes(chromosome):
    """
    Sort chromosomes with a numerical suffix. For example, 'chr1' < 'chr10'.
    
    Args:
    - chromosome (str): Chromosome name as a string.
    
    Returns:
    - tuple: A tuple with the prefix and number part of the chromosome name.
    """
    match = re.match(r'(\D+)(\d+)?', chromosome)
    prefix = match.group(1) if match else ''
    number = int(match.group(2)) if match and match.group(2) else float('inf')
    return (prefix, number)

def read_and_sort_gene_ranges(csv_file_path):
    """
    Reads the CSV file with gene data, expands the gene ranges, and sorts them by chromosome and start position.
    
    Args:
    - csv_file_path (str): The path to the CSV file containing gene data.
    
    Returns:
    - list: A list of dictionaries with expanded gene ranges.
    """
    print("Reading and expanding gene ranges from the CSV file...")
    gene_ranges = []
    try:
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                # Validate necessary fields
                if not all(key in row for key in ["Chromosome", "Start", "End", "GeneName"]):
                    print(f"Skipping invalid row: {row}")
                    continue
                
                try:
                    start = int(row["Start"])
                    end = int(row["End"])
                    gene_length = end - start + 1
                    flank_length = int(4.5 * gene_length)
                    expanded_start = max(0, start - flank_length)  # Ensure no negative start positions
                    expanded_end = end + flank_length

                    gene_range = {
                        "chromosome": row["Chromosome"],
                        "start": expanded_start,
                        "end": expanded_end,
                        "gene_start": start,
                        "gene_end": end,
                        "gene_name": row["GeneName"]
                    }
                    gene_ranges.append(gene_range)
                except ValueError as ve:
                    print(f"Skipping row due to value error: {row} -> {ve}")
            
            # Sort by chromosome and start position
            gene_ranges.sort(key=lambda x: (sort_chromosomes(x["chromosome"]), x["start"]))

        # Print each gene's details
        for idx, gene in enumerate(gene_ranges, start=1):
            print(f"Gene {idx}: {gene['gene_name']}, Chromosome: {gene['chromosome']}, "
                  f"Expanded Range: ({gene['start']}, {gene['end']}), "
                  f"Original Range: ({gene['gene_start']}, {gene['gene_end']})")
        
        print(f"Successfully processed {len(gene_ranges)} gene ranges with expanded boundaries.")
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file_path}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return gene_ranges

def filter_vcf_by_ranges(vcf_file_path, gene_ranges, lines_to_skip=2595):
    """
    Filters a VCF file based on gene ranges and generates a list of matching variants for output.
    
    Args:
    - vcf_file_path (str): Path to the VCF file.
    - gene_ranges (list): List of dictionaries with expanded gene ranges.
    - lines_to_skip (int): Number of lines to skip from the start of the VCF file (for header lines).
    
    Returns:
    - tuple: A tuple containing:
        - header_lines (list): List of header lines to write to the output VCF.
        - csv_rows (list): List of rows to write to the CSV output.
        - vcf_lines (list): List of VCF lines to write to the VCF output.
        - unmatched_genes (dict): A dictionary of genes that did not have any variants matched.
    """
    print("Filtering VCF file based on ranges...")
    csv_rows = []
    vcf_lines = []
    header_lines = []
    total_variants_matched = 0

    # Track genes with no matches
    unmatched_genes = {gene["gene_name"]: True for gene in gene_ranges}

    try:
        with open(vcf_file_path, 'r') as vcf_file:
            # Skip the first lines (header)
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

                # Match against all overlapping gene ranges
                for gene in gene_ranges:
                    if chrom == gene["chromosome"] and gene["start"] <= pos <= gene["end"]:
                        label = gene["gene_name"]
                        if pos < gene["gene_start"] or pos > gene["gene_end"]:
                            label += "_flank"  # Add tag for flanking regions
                        
                        print(f"Match found: Gene={label}, Location={chrom}:{pos}")
                        total_variants_matched += 1

                        # Mark gene as matched
                        unmatched_genes[gene["gene_name"]] = False

                        # Add data for CSV output
                        csv_row = {
                            "GENE": label,
                            "CHROM": vcf_data["CHROM"],
                            "POS": vcf_data["POS"],
                            "ID": vcf_data["ID"],
                            "REF": vcf_data["REF"],
                            "ALT": vcf_data["ALT"],
                            "QUAL": vcf_data["QUAL"],
                            "FILTER": vcf_data["FILTER"],
                            "INFO": vcf_data["INFO"],
                            "FORMAT": vcf_data["FORMAT"],
                            "ORIGINAL_START": gene["gene_start"],
                            "ORIGINAL_END": gene["gene_end"],
                            "EXTENDED_START": gene["start"],
                            "EXTENDED_END": gene["end"]
                        }

                        for idx, sample in enumerate(output_headers[14:]):  # Add genotypes
                            csv_row[sample] = vcf_data["GENOTYPES"][idx] if idx < len(vcf_data["GENOTYPES"]) else "."

                        # Append the CSV row
                        csv_rows.append(csv_row)

                        # Add data for VCF output
                        vcf_lines.append(line.strip())

    except FileNotFoundError:
        print(f"VCF file not found at {vcf_file_path}")
    except Exception as e:
        print(f"Error reading VCF file: {e}")

    print(f"Total variants matched: {total_variants_matched}")
    
    # Print unmatched genes
    print("Genes with no variation found:")
    for gene_name, matched in unmatched_genes.items():
        if matched:
            print(f" - {gene_name}")

    return header_lines, csv_rows, vcf_lines, unmatched_genes

def write_csv(csv_rows):
    """
    Writes the filtered variants to a CSV file.
    
    Args:
    - csv_rows (list): List of rows to write to the CSV file.
    """
    print("Writing the filtered variants to the output CSV file...")
    try:
        with open(output_csv_file_path, 'w', newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=output_headers)
            writer.writeheader()
            writer.writerows(csv_rows)

        print(f"Filtered variants successfully written to {output_csv_file_path}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

def write_vcf(header_lines, vcf_lines):
    """
    Writes the filtered variants to a VCF file.
    
    Args:
    - header_lines (list): List of header lines to write to the VCF file.
    - vcf_lines (list): List of VCF lines to write to the VCF file.
    """
    print("Writing the filtered variants to the output VCF file...")
    try:
        with open(output_vcf_file_path, 'w') as output_file:
            output_file.writelines(header_lines)
            output_file.write("\n".join(vcf_lines) + "\n")
        print(f"Filtered variants successfully written to {output_vcf_file_path}")
    except Exception as e:
        print(f"Error writing VCF file: {e}")

def main():
    """
    Main workflow to filter variants based on gene ranges and output to a chosen format (CSV or VCF).
    """
    print("Starting the filtering workflow...")

    # Step 1: Choose output format
    output_format = input("Enter output format (csv/vcf): ").strip().lower()
    if output_format not in {"csv", "vcf"}:
        print("Invalid output format. Exiting.")
        return

    # Step 2: Read and sort gene ranges
    gene_ranges = read_and_sort_gene_ranges(csv_file_path)
    if not gene_ranges:
        print("No gene ranges available to process. Exiting.")
        return

    # Step 3: Filter the VCF file (first pass)
    header_lines, csv_rows, vcf_lines, unmatched_genes = filter_vcf_by_ranges(vcf_file_path, gene_ranges)

    # Step 4: Re-scan for unmatched genes
    remaining_genes = [gene for gene in gene_ranges if unmatched_genes[gene["gene_name"]]]
    if remaining_genes:
        print(f"Re-scanning for {len(remaining_genes)} unmatched genes...")
        header_lines, csv_rows_second_pass, vcf_lines_second_pass, unmatched_genes_second_pass = filter_vcf_by_ranges(
            vcf_file_path, remaining_genes
        )
        csv_rows.extend(csv_rows_second_pass)
        vcf_lines.extend(vcf_lines_second_pass)

    # Step 5: Write output
    if output_format == "csv":
        write_csv(csv_rows)
    elif output_format == "vcf":
        write_vcf(header_lines, vcf_lines)

    print("Workflow completed.")


if __name__ == "__main__":
    main()

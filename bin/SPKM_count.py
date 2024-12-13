import os
import subprocess
import csv

# Define paths and file names
data_dir = "data/animals"
csv_file = "data/other_animals_RNA.csv"
output_csv = "results/SPKM_counts.csv"
output_header = ["Gene", "hg38_spkm", "galGal6_spkm", "danRer11_spkm", "rn7_spkm", "mm39_spkm", "bosTau9_spkm"]

genomes = ["galGal6", "danRer11", "rn7", "mm39", "bosTau9", "hg38"]

animal_mapping = {
    "galGal6": "Chicken",
    "danRer11": "Zebrafish",
    "rn7": "Rat",
    "mm39": "Mouse",
    "bosTau9": "Cow",
    "hg38": "Human"
}

def calculate_genome_length(genome_path):
    """Calculate genome length using esl-seqstat"""
    try:
        print(f"Calculating genome length for {genome_path}...")
        result = subprocess.run(
            ["esl-seqstat", genome_path],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.split("\n"):
            if line.lstrip().startswith("Total # residues:"):
                total_residues = int(line.split(":")[1].strip().replace(",", "").replace(" ", ""))
                print(f"Genome length for {genome_path}: {total_residues}")
                return total_residues
    except subprocess.CalledProcessError as e:
        print(f"Error calculating genome length for {genome_path}: {e}")
    except ValueError as e:
        print(f"Error parsing genome length for {genome_path}: {e}")
    return None

def calculate_SNP_count(variation_path):
    """Count the number of SNPs in the VCF file"""
    try:
        print(f"Calculating SNP count for {variation_path}...")
        with open(variation_path, "r") as f:
            count = sum(1 for line in f if not line.startswith("#"))
        print(f"SNP count for {variation_path}: {count}")
        return count
    except FileNotFoundError as e:
        print(f"Error opening variation file {variation_path}: {e}")
    return None

def read_and_sort_gene_ranges(csv_file_path, genome):
    """Reads and sorts gene ranges from the provided CSV file"""
    print(f"Reading gene ranges for {genome} from {csv_file_path}...")
    animal = animal_mapping.get(genome, "Unknown")
    gene_ranges = []

    try:
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                loci = row.get(f"{animal}_locus")  # Get the locus for the specific animal
                gene_name = row.get("GeneName", "Unknown")
                if not loci:
                    continue
                if ":" not in loci:
                    continue

                try:
                    chrom, positions = loci.split(":")
                    start, end = map(int, positions.split("-"))

                    # Ensure start is less than end, swap if necessary
                    if start > end:
                        start, end = end, start

                    gene_ranges.append({
                        "chromosome": chrom,
                        "start": start,
                        "end": end,
                        "gene_name": gene_name
                    })
                except ValueError:
                    continue

        gene_ranges.sort(key=lambda x: (x["chromosome"], x["start"]))
        print(f"Successfully sorted {len(gene_ranges)} gene ranges for {genome}.")
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file_path}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return gene_ranges


def filter_vcf_by_ranges(vcf_file_path, gene_ranges):
    """Filters VCF file variants that fall within specified gene ranges"""
    print("Filtering VCF file based on ranges...")
    csv_rows = []
    vcf_lines = []
    header_lines = []
    current_range_idx = 0
    total_variants_matched = 0
    num_variants_matched = {}

    try:
        with open(vcf_file_path, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith("#"):
                    header_lines.append(line)
                    continue

                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = int(fields[1])

                while current_range_idx < len(gene_ranges):
                    gene = gene_ranges[current_range_idx]
                    if chrom < gene["chromosome"] or (chrom == gene["chromosome"] and pos < gene["start"]):
                        break
                    if chrom > gene["chromosome"] or (chrom == gene["chromosome"] and pos > gene["end"]):
                        current_range_idx += 1
                        continue
                    if chrom == gene["chromosome"] and gene["start"] <= pos <= gene["end"]:
                        total_variants_matched += 1
                        if gene["gene_name"] not in num_variants_matched:
                            num_variants_matched[gene["gene_name"]] = 0
                        num_variants_matched[gene["gene_name"]] += 1
                        csv_row = {
                            "GENE": gene["gene_name"],
                            "CHROM": chrom,
                            "POS": pos,
                            "ID": fields[2],
                            "REF": fields[3],
                            "ALT": fields[4],
                            "QUAL": fields[5],
                            "FILTER": fields[6],
                            "INFO": fields[7],
                            "FORMAT": fields[8],
                        }
                        csv_rows.append(csv_row)
                        vcf_lines.append(line.strip())
                        break
    except FileNotFoundError:
        print(f"VCF file not found at {vcf_file_path}")
    except Exception as e:
        print(f"Error reading VCF file: {e}")

    print(f"Total variants matched: {total_variants_matched}")
    return header_lines, csv_rows, vcf_lines, num_variants_matched

def write_csv(csv_rows):
    """Write the filtered variants to CSV file"""
    try:
        with open(output_csv, 'w', newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=["GENE", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"])
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Filtered variants written to {output_csv}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

# Main execution
def main():
    # Initialize results dictionary from the input CSV
    results = {}
    try:
        with open(csv_file, 'r') as csv_file_obj:
            reader = csv.DictReader(csv_file_obj)
            for row in reader:
                gene_name = row.get("GeneName", "Unknown")
                results[gene_name] = {"Gene": gene_name}
    except FileNotFoundError:
        print(f"Input CSV file not found at {csv_file}")
        exit(1)

    for genome in genomes:
        print(f"Processing genome: {genome}")
        genome_path = os.path.join(data_dir, f"{genome}.fa")
        genome_length = calculate_genome_length(genome_path)

        if genome_length is None:
            print(f"Skipping genome {genome} due to genome length calculation error.")
            continue

        variation_path = os.path.join(data_dir, f"{genome}.vcf")
        total_SNP_count = calculate_SNP_count(variation_path)
        print(f"Total SNP count for {genome}: {total_SNP_count}")

        if total_SNP_count is None:
            print(f"Skipping genome {genome} due to SNP count calculation error.")
            continue

        gene_ranges = read_and_sort_gene_ranges(csv_file, genome)
        header_lines, csv_rows, vcf_lines, num_variants_matched = filter_vcf_by_ranges(variation_path, gene_ranges)

        # Calculate SPKM for each gene
        for gene_name, SNP_count in num_variants_matched.items():
            # Ensure that SNP_count is an integer, not a dictionary
            if isinstance(SNP_count, int):  # Check if it's an integer
                SPKM = (10**6) * (1000 * SNP_count / genome_length) / total_SNP_count
                results[gene_name][f"{genome}_spkm"] = SPKM
            else:
                print(f"Skipping {gene_name} due to invalid SNP_count: {SNP_count}")

    # Write results to CSV
    write_csv(csv_rows)

if __name__ == "__main__":
    main()

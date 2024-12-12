import os
import subprocess
import csv
import pandas as pd

# Define paths and file names
data_dir = "data/animals"
csv_file = "data/other_animals_RNA.csv"
output_csv = "results/SPKM_counts.csv"

output_header = ["Gene", "hg38_spkm", "galGal6_spkm", "danRer11_spkm", "rn7_spkm", "mm39_spkm", "bosTau9_spkm"]

genomes = ["galGal6", "danRer11", "rn7", "mm39", "bosTau9", "hg38"]

def calculate_genome_length(genome_path):
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
    print(f"Reading gene ranges for {genome} from {csv_file_path}...")
    animal_mapping = {
        "galGal6": "Chicken",
        "danRer11": "Zebrafish",
        "rn7": "Rat",
        "mm39": "Mouse",
        "bosTau9": "Cow",
        "hg38": "Human"
    }
    animal = animal_mapping.get(genome, "Unknown")
    gene_ranges = []

    try:
        with open(csv_file_path, 'r') as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                loci = row.get(f"{animal}_locus")
                gene_name = row.get("GeneName", "Unknown")
                if not loci:
                    continue
                if ":" not in loci:
                    continue

                try:
                    chrom, positions = loci.split(":")
                    start, end = map(int, positions.split("-"))
                    if start < end:
                        gene_ranges.append({
                            "chromosome": chrom,
                            "start": start,
                            "end": end,
                            "gene_name": gene_name
                        })
                except ValueError:
                    continue

        gene_ranges.sort(key=lambda x: (x["chromosome"], x["start"]))
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file_path}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return gene_ranges

def filter_vcf_by_ranges(vcf_file_path, gene_ranges):
    current_range_idx = 0
    num_variants_matched = {}

    try:
        with open(vcf_file_path, 'r') as vcf_file:
            for line in vcf_file:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                chrom = fields[0].lstrip("chr")
                pos = int(fields[1])

                while current_range_idx < len(gene_ranges):
                    gene = gene_ranges[current_range_idx]
                    gene_chrom = gene["chromosome"].lstrip("chr")

                    if chrom < gene_chrom or (chrom == gene_chrom and pos < gene["start"]):
                        break

                    if chrom > gene_chrom or (chrom == gene_chrom and pos > gene["end"]):
                        current_range_idx += 1
                        continue

                    if chrom == gene_chrom and gene["start"] <= pos <= gene["end"]:
                        gene_name = gene["gene_name"]
                        num_variants_matched[gene_name] = num_variants_matched.get(gene_name, 0) + 1
                        break
        
    except FileNotFoundError:
        print(f"VCF file not found at {vcf_file_path}")
    except Exception as e:
        print(f"Error reading VCF file: {e}")
    return num_variants_matched

# Initialize results dictionary with all genes from the input CSV
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
    num_variants_matched = filter_vcf_by_ranges(variation_path, gene_ranges)

    for gene_name, SNP_count in num_variants_matched.items():
        SPKM = (10**6) * (1000 * SNP_count / genome_length) / total_SNP_count
        results[gene_name][f"{genome}_spkm"] = SPKM

# Convert results to DataFrame
print("Converting results to DataFrame...")
df = pd.DataFrame.from_dict(results, orient="index").reindex(columns=output_header, fill_value="NA")
df.to_csv(output_csv, index=False)
print(f"Processing complete. Results saved to {output_csv}.")

import os
import subprocess
import csv

# Define paths and file names
data_dir = "data/animals"
csv_file = "data/animals_var.csv"
output_csv = "results/SPKM_counts.csv"
output_header = ["Gene", "hg38_spkm", "galGal6_spkm", "mm39_spkm", "rn7_spkm", "danRer11_spkm", "bosTau9_spkm"]

genomes = ["hg38", "galGal6", "mm39", "rn7", "danRer11", "bosTau9"]

# Map genomes to their common names
animal_mapping = {
    "galGal6": "Chicken",
    "danRer11": "Zebrafish",
    "rn7": "Rat",
    "mm39": "Mouse",
    "bosTau9": "Cow",
    "hg38": "Human"
}

def get_genome_data():
    """Precompute genome lengths and total SNP counts for each genome"""
    genome_data = {}
    for genome in genomes:
        genome_path = os.path.join(data_dir, f"{genome}.fa")
        genome_length = calculate_genome_length(genome_path)
        total_snp_count = calculate_total_snp_count(genome)
        
        genome_data[genome] = {
            "length": genome_length,
            "total_snp_count": total_snp_count
        }
    return genome_data

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

def calculate_total_snp_count(genome):
    """Calculate the total variation count for each genome from input CSV"""
    total_snp_count = 0
    try:
        with open(csv_file, 'r') as csv_file_obj:
            reader = csv.DictReader(csv_file_obj)
            for row in reader:
                total_snp_count += int(row.get(f"{animal_mapping[genome]}_var", 0))
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return total_snp_count

def read_gene_variation_counts():
    """Reads gene variation counts from the provided CSV file"""
    print(f"Reading gene variation counts from {csv_file}...")
    gene_variations = {}

    try:
        with open(csv_file, 'r') as csv_file_obj:
            reader = csv.DictReader(csv_file_obj)
            for row in reader:
                gene_name = row.get("Gene", "Unknown")
                if gene_name != "Unknown":
                    # Read variation counts for each genome
                    gene_variations[gene_name] = {
                        "hg38": int(row.get("Human_var", 0)),
                        "galGal6": int(row.get("Chicken_var", 0)),
                        "danRer11": int(row.get("Zebrafish_var", 0)),
                        "rn7": int(row.get("Rat_var", 0)),
                        "mm39": int(row.get("Mouse_var", 0)),
                        "bosTau9": int(row.get("Cow_var", 0)),
                    }
    except FileNotFoundError:
        print(f"CSV file not found at {csv_file}")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
    return gene_variations

def write_results_to_csv(results):
    """Write the calculated SPKM results to CSV file"""
    try:
        with open(output_csv, 'w', newline='') as output_file:
            writer = csv.DictWriter(output_file, fieldnames=output_header)
            writer.writeheader()
            writer.writerows(results)
        print(f"SPKM results written to {output_csv}")
    except Exception as e:
        print(f"Error writing CSV file: {e}")

def main():
    # Precompute genome lengths and total SNP counts for each genome
    genome_data = get_genome_data()

    # Read variation counts from input CSV
    gene_variations = read_gene_variation_counts()

    # Initialize results list
    results = []

    # Process each gene and calculate SPKM for each genome
    for gene_name, variations in gene_variations.items():
        gene_result = {"Gene": gene_name}

        # Process each genome and calculate SPKM
        for genome in genomes:
            genome_data_for_genome = genome_data.get(genome)
            if not genome_data_for_genome or genome_data_for_genome["length"] is None or genome_data_for_genome["total_snp_count"] is None:
                print(f"Skipping genome {genome} due to missing data.")
                gene_result[f"{genome}_spkm"] = 0
                continue

            genome_length = genome_data_for_genome["length"]
            total_snp_count = genome_data_for_genome["total_snp_count"]
            variation_count = variations.get(genome, 0)

            # Calculate SPKM for this gene in the current genome
            if total_snp_count > 0 and genome_length > 0:
                SPKM = (10**6) * (1000 * variation_count / genome_length) / total_snp_count
                gene_result[f"{genome}_spkm"] = SPKM
            else:
                gene_result[f"{genome}_spkm"] = 0  # Set to 0 if any value is invalid or 0

        results.append(gene_result)

    # Write the results to a CSV
    write_results_to_csv(results)

if __name__ == "__main__":
    main()

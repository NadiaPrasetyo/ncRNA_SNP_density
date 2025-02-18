import csv
import os
import pandas as pd

# Input filenames
csv_file = "data/SNP-densities-and-RNA.csv"
sam_file = "data/aligned_genome_reads.sam"
output_file = "results/adjusted_alignments.csv"

# Read reference CSV file into a dictionary
reference_data = {}
print("Reading reference CSV file...")
with open(csv_file, newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_name = row["GeneName"].strip()
        reference_data[gene_name] = {
            "chromosome": row["Chromosome"],
            "start": int(row["Start"])
        }
print(f"Loaded {len(reference_data)} gene references.")

# Process SAM file and adjust alignments
output_data = []
missing_genes = set()
total_lines = sum(1 for _ in open(sam_file, "r", encoding="utf-8"))
print(f"Processing SAM file ({total_lines} lines)...")
reference_locations = {}  # Store reference genome locations for adjustment

with open(sam_file, "r", encoding="utf-8") as samfile:
    for line_number, line in enumerate(samfile, start=1):
        if line_number % 1000 == 0:
            print(f"Processed {line_number}/{total_lines} lines...")
        
        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue  # Skip malformed lines
        
        filename_parts = fields[0].replace(".sam", "").rsplit("-", 1)  # Split only at the last hyphen
        if len(filename_parts) < 2:
            continue  # Skip malformed filenames
        
        gene_name, genome = filename_parts[0], filename_parts[1]
        if gene_name not in reference_data:
            missing_genes.add(gene_name)
            continue  # Skip genes not in reference CSV
        
        chromosome = reference_data[gene_name]["chromosome"]
        reference_start = reference_data[gene_name]["start"]
        length = len(fields[10])
        alignment_score = int(fields[5])
        
        genome_start = int(fields[4])
        
        if genome == "HG00438":
            reference_locations[gene_name] = genome_start
            start_location = reference_start
        else:
            if gene_name in reference_locations:
                hg00438_start = reference_locations[gene_name]
                start_location = reference_start + (genome_start - hg00438_start)
            else:
                start_location = genome_start  # Fallback if HG00438 is missing
        
        output_data.append([gene_name, genome, chromosome, start_location, length, alignment_score])

print(f"Finished processing SAM file. Missing genes in some genomes: {missing_genes}")

# Save results to CSV
output_df = pd.DataFrame(output_data, columns=["Gene", "Genome", "Chromosome", "Start Location", "Length", "Alignment Score"])
output_df.to_csv(output_file, index=False)

print(f"Processed SAM file and saved results to {output_file}")

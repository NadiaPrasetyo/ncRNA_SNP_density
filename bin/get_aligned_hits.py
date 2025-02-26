import csv
import re
from collections import defaultdict

def parse_sam_file(sam_file):
    gene_hits = defaultdict(lambda: defaultdict(int))
    
    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith('@'):  # Skip header lines
                continue
            
            fields = line.strip().split('\t')
            filename = fields[0]  # Extract filename (gene-genome.sam)
            
            match = re.match(r"(.+)-(.+).sam", filename)
            if match:
                gene, genome = match.groups()
                gene_hits[gene][genome] += 1
            else:
                print(f"Warning: Unable to parse filename {filename}")
    
    return gene_hits

def write_csv(output_file, gene_hits):
    genomes = sorted({genome for hits in gene_hits.values() for genome in hits})
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        header = ['Gene name'] + genomes
        writer.writerow(header)
        
        for gene, hits in sorted(gene_hits.items()):
            row = [gene] + [hits.get(genome, 0) for genome in genomes]
            writer.writerow(row)

def main():
    sam_file = "data/aligned_genome_reads.sam"
    output_file = "results/gene_hits_checkboard.csv"
    
    gene_hits = parse_sam_file(sam_file)
    write_csv(output_file, gene_hits)
    
    print(f"CSV file '{output_file}' has been created successfully.")

if __name__ == "__main__":
    main()

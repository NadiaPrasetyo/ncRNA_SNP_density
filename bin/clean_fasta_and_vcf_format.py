import csv
import re

def load_gene_mapping(csv_file):
    gene_mapping = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_mapping.append({
                'CHROM': row['Chromosome'],
                'START': int(row['Start']),
                'END': int(row['End']),
                'GENE': row['GeneName']
            })
    return gene_mapping

def map_gene(chrom, start, end, gene_mapping):
    for entry in gene_mapping:
        if entry['CHROM'] == chrom and entry['START'] <= start <= entry['END']:
            return entry['GENE']
        if entry['CHROM'] == chrom and entry['START'] <= end <= entry['END']:
            return entry['GENE']
    return "Unknown"

def reformat_fasta(input_fasta, output_fasta, gene_mapping):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract chromosome and position range
                header_match = re.search(r'NC_(\d+)\.\d+:(\d+)-(\d+)', line)
                if header_match:
                    chrom = f"chr{int(header_match.group(1))}"  # Convert NC_000016.10 to chr16
                    start = int(header_match.group(2))
                    end = int(header_match.group(3))
                    gene = map_gene(chrom, start, end, gene_mapping)
                    outfile.write(f">{gene}\n")
                else:
                    outfile.write(f">Unknown\n")
            else:
                outfile.write(line)

def reformat_vcf_with_gene(input_vcf, output_vcf, gene_mapping):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)
            elif line.startswith('#CHROM'):
                # Add the GENE column to the header
                outfile.write("#GENE\t" + line)
            else:
                fields = line.strip().split('\t')
                chrom, pos = fields[0], int(fields[1])
                gene = map_gene(chrom, pos, pos, gene_mapping)  # Pass start and end as the same value for a single position
                outfile.write(f"{gene}\t" + line)


# File paths
input_fasta = "data/combined_sequences.fa"
output_fasta = "data/CLEAN_combined_sequences.fa"
input_vcf = "data/filtered_variants.vcf"
output_vcf = "data/CLEAN_filtered_variants.vcf"
gene_csv = "data/SNP-densities-and-RNA.csv"

# Load gene mapping
gene_mapping = load_gene_mapping(gene_csv)

# Run the reformatting functions
# Run the FASTA reformatting
reformat_fasta(input_fasta, output_fasta, gene_mapping)
reformat_vcf_with_gene(input_vcf, output_vcf, gene_mapping)

print("Reformatting complete.")

import os
import re
import csv

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

def reformat_fasta(input_fasta, gene_mapping):
    # Read the content of the fasta file and write it back after modifications
    with open(input_fasta, 'r') as infile:
        lines = infile.readlines()
    
    with open(input_fasta, 'w') as outfile:
        for line in lines:
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
    """
    Reformat the VCF file to use gene names as #CHROM and relative positions as POS.
    """
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith('##'):
                # Preserve metadata lines
                outfile.write(line)
            elif line.startswith('#CHROM'):
                # Add GENE as the #CHROM header
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                chrom, pos = fields[0], int(fields[1])
                gene, gene_start = None, None

                # Identify the gene and calculate relative position
                for entry in gene_mapping:
                    if entry['CHROM'] == chrom and entry['START'] <= pos <= entry['END']:
                        gene = entry['GENE']
                        gene_start = entry['START']
                        break

                if gene:
                    relative_pos = pos - gene_start + 1  # Compute 1-based relative position
                    fields[0] = gene  # Replace CHROM with gene name
                    fields[1] = str(relative_pos)  # Update POS to relative position
                    outfile.write("\t".join(fields) + "\n")
                else:
                    # If no gene found, skip the variant
                    print(f"Warning: No gene found for variant {chrom}:{pos}")


# File paths
input_vcf = "data/filtered_variants.vcf"
output_vcf = "data/CLEAN_filtered_variants.vcf"
gene_csv = "data/SNP-densities-and-RNA.csv"

# Load gene mapping
gene_mapping = load_gene_mapping(gene_csv)

# Process all .fasta files in the 'FASTA' directory and rename them to .fa
fasta_directory = "data/FASTA"
for filename in os.listdir(fasta_directory):
    if filename.endswith(".fasta"):
        input_fasta = os.path.join(fasta_directory, filename)
        
        # Reformat the FASTA file
        reformat_fasta(input_fasta, gene_mapping)
        
        # Rename the file to .fa
        new_filename = filename.replace(".fasta", ".fa")
        new_filepath = os.path.join(fasta_directory, new_filename)
        
        # Rename the file
        os.rename(input_fasta, new_filepath)
        print(f"Renamed {filename} to {new_filename}")

# Run the VCF reformatting (you can modify the VCF as well if needed)
reformat_vcf_with_gene(input_vcf, output_vcf, gene_mapping)

print("Reformatting and renaming complete.")
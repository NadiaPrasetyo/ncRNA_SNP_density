import csv
import re
from collections import defaultdict

def parse_csv(csv_file):
    """Parse the SNP-densities-and-RNA.csv file to extract gene position data."""
    gene_data = []
    with open(csv_file, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_data.append({
                'Chromosome': row['Chromosome'],
                'Start': int(row['Start']),
                'End': int(row['End']),
                'GeneName': row['GeneName'].strip()  # Remove extra whitespace
            })
    return gene_data

def parse_sam(sam_file):
    """Parse the SAM file to extract mapping information."""
    sam_data = []
    with open(sam_file, mode='r') as file:
        for line in file:
            if line.startswith('@'):  # Skip headers
                continue
            parts = line.split('\t')
            match = re.search(r'chr(\d+)', parts[0])  # Extract chromosome info from filename
            gene = parts[1].strip()  # Remove extra whitespace
            mapped_chromosome = f"chr{match.group(1)}" if match else None
            mapped_position = int(parts[4]) if parts[4].isdigit() else None
            
            sam_data.append({
                'Gene': gene,
                'Mapped Chromosome': mapped_chromosome,
                'Mapped Position': mapped_position,
                'CIGAR': parts[6],
            })
    return sam_data

def calculate_length(cigar):
    """Calculate the alignment length from the CIGAR string."""
    matches = re.findall(r'(\d+)[MIDNSHP=X]', cigar)
    return sum(int(m) for m in matches)

def summarize_data(sam_data, gene_data):
    """Summarize data based on gene names."""
    summary = []
    for sam_entry in sam_data:
        gene_name = sam_entry['Gene']
               
        matched_gene = next((gene for gene in gene_data if gene['GeneName'] == gene_name), None)
        if matched_gene:
            summary.append({
                'Gene': matched_gene['GeneName'],
                'Chromosome': matched_gene['Chromosome'],
                'Position': f"{matched_gene['Start']}-{matched_gene['End']}",
                'Mapped Chromosome': sam_entry['Mapped Chromosome'],
                'Mapped Position': sam_entry['Mapped Position'],
                'Length': calculate_length(sam_entry['CIGAR']),
            })
        else:
            print(f"No match for gene {gene_name}")
    return summary

def concise_summary(sam_data, gene_data):
    """Create a concise summary with Gene, Chromosome, and a list of Mapped Chromosomes."""
    summary = defaultdict(set)  # Use a dictionary to collect mapped chromosomes for each gene
    
    for sam_entry in sam_data:
        gene_name = sam_entry['Gene']
        matched_gene = next((gene for gene in gene_data if gene['GeneName'] == gene_name), None)
        if matched_gene:
            summary[gene_name].add(sam_entry['Mapped Chromosome'])
        else:
            print(f"No match for gene {gene_name}")
    
    concise_list = [
        {
            'Gene': gene,
            'Chromosome': next((gene_data_item['Chromosome'] for gene_data_item in gene_data if gene_data_item['GeneName'] == gene), "N/A"),
            'Mapped Chromosomes': ", ".join(sorted(mapped_chromosomes))  # Sort and join mapped chromosomes
        }
        for gene, mapped_chromosomes in summary.items()
    ]
    return concise_list

def save_summary(summary, output_file):
    """Save the summarized data to a CSV file."""
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Gene', 'Chromosome', 'Position', 'Mapped Chromosome', 'Mapped Position', 'Length'])
        writer.writeheader()
        writer.writerows(summary)
        
def save_concise_summary(concise_list, output_file):
    """Save the concise summary to a CSV file."""
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Gene', 'Chromosome', 'Mapped Chromosomes'])
        writer.writeheader()
        writer.writerows(concise_list)


# File paths
csv_file = 'data/SNP-densities-and-RNA.csv'
sam_file = 'data/aligned_reads.sam'
output_file = 'results/aligned_genes_summary.csv'
concise_summary_file = 'results/concise_aligned_genes_summary.csv'

# Processing
gene_data = parse_csv(csv_file)
sam_data = parse_sam(sam_file)
summary = summarize_data(sam_data, gene_data)
concise_list = concise_summary(sam_data, gene_data)
save_summary(summary, output_file)
save_concise_summary(concise_list, concise_summary_file)

print(f"Summary saved to {output_file}")

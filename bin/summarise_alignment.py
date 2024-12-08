import csv
import re
from collections import defaultdict

def parse_csv(csv_file):
    """
    Parse the SNP-densities-and-RNA.csv file to extract gene position data.
    
    Args:
        csv_file (str): Path to the CSV file containing gene position data.
    
    Returns:
        list: A list of dictionaries containing gene data with Chromosome, Start, End, and GeneName.
    """
    gene_data = []
    try:
        with open(csv_file, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                try:
                    gene_data.append({
                        'Chromosome': row['Chromosome'],
                        'Start': int(row['Start']),
                        'End': int(row['End']),
                        'GeneName': row['GeneName'].strip()  # Remove extra whitespace
                    })
                except (KeyError, ValueError) as e:
                    print(f"Error processing row in CSV file: {e}")
    except FileNotFoundError:
        print(f"Error: The file '{csv_file}' was not found.")
    except IOError as e:
        print(f"Error reading file '{csv_file}': {e}")
    return gene_data

def parse_sam(sam_file):
    """
    Parse the SAM file to extract mapping information.
    
    Args:
        sam_file (str): Path to the SAM file containing mapping data.
    
    Returns:
        list: A list of dictionaries containing gene mapping data.
    """
    sam_data = []
    try:
        with open(sam_file, mode='r') as file:
            for line in file:
                if line.startswith('@'):  # Skip headers
                    continue
                try:
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
                except IndexError as e:
                    print(f"Error parsing line in SAM file: {e}")
    except FileNotFoundError:
        print(f"Error: The file '{sam_file}' was not found.")
    except IOError as e:
        print(f"Error reading file '{sam_file}': {e}")
    return sam_data

def calculate_length(cigar):
    """
    Calculate the alignment length from the CIGAR string.
    
    Args:
        cigar (str): CIGAR string from the SAM file.
    
    Returns:
        int: Total alignment length.
    """
    try:
        matches = re.findall(r'(\d+)[MIDNSHP=X]', cigar)
        return sum(int(m) for m in matches)
    except Exception as e:
        print(f"Error calculating length from CIGAR: {e}")
        return 0

def summarize_data(sam_data, gene_data):
    """
    Summarize data by combining SAM and gene information.
    
    Args:
        sam_data (list): List of SAM data dictionaries.
        gene_data (list): List of gene data dictionaries.
    
    Returns:
        list: A summary of aligned genes.
    """
    summary = []
    for sam_entry in sam_data:
        try:
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
        except Exception as e:
            print(f"Error processing SAM entry {sam_entry}: {e}")
    return summary

def concise_summary(sam_data, gene_data):
    """
    Create a concise summary of aligned genes.
    
    Args:
        sam_data (list): List of SAM data dictionaries.
        gene_data (list): List of gene data dictionaries.
    
    Returns:
        list: A concise summary with gene and mapped chromosomes.
    """
    summary = defaultdict(set)
    for sam_entry in sam_data:
        try:
            gene_name = sam_entry['Gene']
            matched_gene = next((gene for gene in gene_data if gene['GeneName'] == gene_name), None)
            if matched_gene:
                summary[gene_name].add(sam_entry['Mapped Chromosome'])
            else:
                print(f"No match for gene {gene_name}")
        except Exception as e:
            print(f"Error processing SAM entry {sam_entry}: {e}")
    
    concise_list = [
        {
            'Gene': gene,
            'Chromosome': next((g['Chromosome'] for g in gene_data if g['GeneName'] == gene), "N/A"),
            'Mapped Chromosomes': ", ".join(sorted(mapped_chromosomes))  # Sort and join mapped chromosomes
        }
        for gene, mapped_chromosomes in summary.items()
    ]
    return concise_list

def save_summary(summary, output_file):
    """
    Save the detailed summary to a CSV file.
    
    Args:
        summary (list): List of summarized gene data.
        output_file (str): Path to save the output CSV file.
    """
    try:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['Gene', 'Chromosome', 'Position', 'Mapped Chromosome', 'Mapped Position', 'Length'])
            writer.writeheader()
            writer.writerows(summary)
    except IOError as e:
        print(f"Error writing file '{output_file}': {e}")

def save_concise_summary(concise_list, output_file):
    """
    Save the concise summary to a CSV file.
    
    Args:
        concise_list (list): List of concise summarized gene data.
        output_file (str): Path to save the output CSV file.
    """
    try:
        with open(output_file, mode='w', newline='') as file:
            writer = csv.DictWriter(file, fieldnames=['Gene', 'Chromosome', 'Mapped Chromosomes'])
            writer.writeheader()
            writer.writerows(concise_list)
    except IOError as e:
        print(f"Error writing file '{output_file}': {e}")


# File paths
csv_file = 'data/SNP-densities-and-RNA.csv'
sam_file = 'data/aligned_reads.sam'
output_file = 'results/aligned_genes_summary.csv'
concise_summary_file = 'results/concise_aligned_genes_summary.csv'

# Processing
try:
    gene_data = parse_csv(csv_file)
    sam_data = parse_sam(sam_file)
    summary = summarize_data(sam_data, gene_data)
    concise_list = concise_summary(sam_data, gene_data)
    save_summary(summary, output_file)
    save_concise_summary(concise_list, concise_summary_file)
    print(f"Summary saved to {output_file}")
    print(f"Concise summary saved to {concise_summary_file}")
except Exception as e:
    print(f"An unexpected error occurred during processing: {e}")

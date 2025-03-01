import pandas as pd
import re
from collections import defaultdict

def parse_snp_data(file_path):
    snp_counts = defaultdict(lambda: defaultdict(int))
    with open(file_path, 'r') as f:
        content = f.read()
    
    entries = content.split('----------------------------------------')
    for entry in entries:
        entry = entry.strip()
        if not entry:
            continue
        
        gene_match = re.search(r'Gene: (.+)', entry)
        ref_match = re.search(r'Reference Allele: (.+)', entry)
        alt_match = re.search(r'Alternate Alleles: (.+)', entry)
        variation_class_match = re.search(r'Variation Class: (.+)', entry)
        
        if gene_match and ref_match and alt_match and variation_class_match:
            gene = gene_match.group(1)
            ref_allele = ref_match.group(1)
            alt_alleles = alt_match.group(1).split(', ')
            variation_class = variation_class_match.group(1)

            # Only process if the variation class is 'snv'
            if variation_class == 'snv':
                for alt in alt_alleles:
                    mutation = f"{ref_allele}->{alt}"
                    if is_single_base_mutation(mutation):
                        snp_counts[gene][mutation] += 1
    
    return snp_counts

def parse_pangenome_data(file_path):
    pangenome_counts = defaultdict(lambda: defaultdict(int))
    df = pd.read_csv(file_path, dtype=str)
    df = df[df['TYPE'] == 'snp']
    
    for _, row in df.iterrows():
        gene = row['GENE']
        mutations = row['SNP_TYPE'].split(',')  # Handle multiple mutations
        
        if pd.notna(row['GENOME_DETAILS']):
            genome_entries = row['GENOME_DETAILS'].split('; ')
            
            for mutation in mutations:
                mutation = mutation.strip().replace('>', '->')  # Convert format
                if is_single_base_mutation(mutation):
                    count = 0
                    
                    for genome in genome_entries:
                        if "Homozygous" in genome:
                            count += 2
                        elif "Heterozygous" in genome:
                            count += 1
                    
                    pangenome_counts[gene][mutation] += count
    
    return pangenome_counts

def is_single_base_mutation(mutation):
    # Split the mutation string into ref and alt alleles
    ref, alt = mutation.split('->')
    
    # Check if both ref and alt alleles consist of a single base
    return len(ref) == 1 and len(alt) == 1

def merge_counts(snp_data_counts, pangenome_counts):
    merged_counts = defaultdict(lambda: defaultdict(lambda: {'SNP155': 0, 'Pangenome': 0}))
    
    for gene, mutations in snp_data_counts.items():
        for mutation, count in mutations.items():
            merged_counts[gene][mutation]['SNP155'] += count
    
    for gene, mutations in pangenome_counts.items():
        for mutation, count in mutations.items():
            merged_counts[gene][mutation]['Pangenome'] += count
    
    return merged_counts

def write_output(merged_counts, output_file):
    output_data = []
    for gene, mutations in merged_counts.items():
        for mutation, counts in mutations.items():
            output_data.append([gene, mutation, counts['SNP155'], counts['Pangenome']])
    
    df_output = pd.DataFrame(output_data, columns=["Gene", "Mutation", "SNP155_Count", "Pangenome_Count"])
    df_output.to_csv(output_file, index=False)

def main():
    snp_data_counts = parse_snp_data("data/SNP_data.txt")
    pangenome_counts = parse_pangenome_data("data/pangenome_summary.csv")
    merged_counts = merge_counts(snp_data_counts, pangenome_counts)
    write_output(merged_counts, "data/snp_frequencies.csv")
    print("Processing complete. Output saved to data/snp_frequencies.csv")

if __name__ == "__main__":
    main()

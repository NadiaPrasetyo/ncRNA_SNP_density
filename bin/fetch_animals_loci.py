import pandas as pd
import requests
import time

# Function to query UCSC API and retrieve locus data
def get_locus_data(gene_name, genome):
    url = f"http://api.genome.ucsc.edu/search?search={gene_name}&genome={genome}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        # Extract the chromosome location
        for item in data.get('results', []):
            if 'position' in item:
                return item['position']
        return ""
    except Exception as e:
        print(f"Error querying {gene_name} for {genome}: {e}")
        return ""

# Function to process the file and add new locus data
def add_locus_data(input_file, output_file):
    # Load input data
    df = pd.read_csv(input_file)

    # Initialize new columns for each species
    df['Human_locus'] = ""
    df['Chicken_locus'] = ""
    df['Mouse_locus'] = ""
    df['Rat_locus'] = ""
    df['Zebrafish_locus'] = ""
    df['Cow_locus'] = ""

    # Define genome identifiers
    genomes = {
        'Human_locus': 'hg38',
        'Chicken_locus': 'galGal6',
        'Mouse_locus': 'mm39',
        'Rat_locus': 'rn7',
        'Zebrafish_locus': 'danRer11',
        'Cow_locus': 'bosTau9'
    }

    # Query UCSC API for each gene and genome
    for index, row in df.iterrows():
        gene_name = row['GeneName']
        
        for column, genome in genomes.items():
            locus = get_locus_data(gene_name, genome)
            df.at[index, column] = locus
            time.sleep(0.5)  # To avoid overwhelming the API

    # Save updated data to the output file
    df.to_csv(output_file, index=False)

# Example usage
input_file = "data/SNP-densities-and-RNA.csv"
output_file = "data/other_animals_RNA.csv"
add_locus_data(input_file, output_file)

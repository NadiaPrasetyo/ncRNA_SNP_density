import pandas as pd
import requests
import time

# Function to query UCSC API and retrieve locus data
def get_locus_data(gene_name, genome):
    url = f"http://api.genome.ucsc.edu/search?search={gene_name}&genome={genome}"
    print(f"Querying {url}")  # Diagnostic print
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        # Extract the chromosome location from the nested matches
        for item in data.get('positionMatches', []):
            for match in item.get('matches', []):
                if 'position' in match:
                    print(f"Found position for {gene_name} in {genome}: {match['position']}")  # Diagnostic print
                    return match['position']
        print(f"No position found for {gene_name} in {genome}")  # Diagnostic print
        return ""
    except Exception as e:
        print(f"Error querying {gene_name} for {genome}: {e}")  # Diagnostic print
        return ""

# Function to process the file and add new locus data
def add_locus_data(input_file, output_file):
    print(f"Loading input file: {input_file}")  # Diagnostic print
    # Load input data
    df = pd.read_csv(input_file)

    # Extract Human locus from existing fields
    df['Human_locus'] = df.apply(lambda row: f"{row['Chromosome']}:{row['Start']}-{row['End']}", axis=1)

    # Initialize new columns for other species
    df['Chicken_locus'] = ""
    df['Mouse_locus'] = ""
    df['Rat_locus'] = ""
    df['Zebrafish_locus'] = ""
    df['Cow_locus'] = ""

    # Define genome identifiers
    genomes = {
        'Chicken_locus': 'galGal6',
        'Mouse_locus': 'mm39',
        'Rat_locus': 'rn7',
        'Zebrafish_locus': 'danRer11',
        'Cow_locus': 'bosTau9'
    }

    # Query UCSC API for each gene and genome
    for index, row in df.iterrows():
        gene_name = row['GeneName']
        print(f"Processing gene {gene_name} (row {index + 1}/{len(df)})")  # Diagnostic print
        
        for column, genome in genomes.items():
            print(f"  Querying {column} for genome {genome}")  # Diagnostic print
            locus = get_locus_data(gene_name, genome)
            df.at[index, column] = locus
            time.sleep(0.5)  # To avoid overwhelming the API

    # Save updated data to the output file
    print(f"Saving updated data to {output_file}")  # Diagnostic print
    df.to_csv(output_file, index=False)
    print("Processing complete.")  # Diagnostic print

# Example usage
input_file = "data/filtered_ncRNA_list.csv"
output_file = "data/other_animals_RNA.csv"
add_locus_data(input_file, output_file)

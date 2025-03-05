import pandas as pd
import requests
import random
import json

# Define the GraphQL endpoint
gnomad_api_url = "https://gnomad.broadinstitute.org/api"

# Function to fetch data from the GraphQL API
def fetch_data_from_gnomad(chrom, start, stop):
    # Remove the 'chr' prefix if present (gnomAD expects just the number)
    chrom = chrom.replace('chr', '')
    
    # Define the GraphQL query
    query = {
        "query": """
        query {
          region(chrom: "%s", start: %d, stop: %d, reference_genome: GRCh38) {
            chrom
            variants(dataset: gnomad_r4) {
              variant_id
              hgvsc
              genome {
                ac
                an
              }
              exome {
                ac
                an
              }
            }
          }
        }
        """ % (chrom, start, stop)
    }
    
    # Debugging: print the query to see the structure
    print(f"Sending query for {chrom}:{start}-{stop}...")
    
    # Send the POST request to the GraphQL API
    response = requests.post(gnomad_api_url, json=query)
    
    # Check if the response is successful
    if response.status_code == 200:
        return response.json()
    else:
        # Print error diagnostics
        print(f"Error fetching data for {chrom}:{start}-{stop} - Status Code: {response.status_code}")
        print(f"Response Body: {response.text}")
        return None

# Read the CSV file
csv_file = 'data/SNP_RNA_GC.csv'
df = pd.read_csv(csv_file)

# Select the first 20 rows
first_20_rows = df.head(20)

# Select 100 random rows starting from line 124 (index 123)
random_rows = df.iloc[123:].sample(n=100)

# Combine the selected rows
selected_rows = pd.concat([first_20_rows, random_rows])

# Prepare a list to hold the result data
results = []

# Loop through the selected rows and fetch data
total_rows = len(selected_rows)
for idx, (_, row) in enumerate(selected_rows.iterrows(), start=1):
    chrom = row['Chromosome']
    start = int(row['Start'])
    stop = int(row['End'])
    gene_name = row['GeneName']
    
    # Print progress
    print(f"Progress: {idx}/{total_rows} - Processing {chrom}:{start}-{stop}")
    
    # Fetch the data from gnomAD
    result = fetch_data_from_gnomad(chrom, start, stop)
    
    # Check if the result is not None and process it
    if result:
        variants = result.get('data', {}).get('region', {}).get('variants', [])
        if not variants:
            print(f"No variants found for {chrom}:{start}-{stop}")
        for variant in variants:
            variant_id = variant.get('variant_id')
            hgvsc = variant.get('hgvsc')
            
            # Determine whether to use genome or exome data
            genome_data = variant.get('genome')
            exome_data = variant.get('exome')
            
            if genome_data:
                ac = genome_data.get('ac')
                an = genome_data.get('an')
            elif exome_data:
                ac = exome_data.get('ac')
                an = exome_data.get('an')
            else:
                continue  # Skip if neither genome nor exome data is available
            
            allele_freq = ac / an if an else None
            
            # Append the data to the results list
            results.append({
                'gene_name': gene_name,
                'chrom': chrom,
                'start': start,
                'stop': stop,
                'variation_id': variant_id,
                'variation_consequence': hgvsc,
                'allele_count': ac,
                'allele_num': an,
                'allele_freq': allele_freq
            })

# Create a DataFrame from the results and save it to a CSV file
output_df = pd.DataFrame(results)
output_csv_file = 'data/gnomad_gene_rand_freq.csv'
output_df.to_csv(output_csv_file, index=False)

print(f"Data has been saved to {output_csv_file}")

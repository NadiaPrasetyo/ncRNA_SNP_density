import pandas as pd
import requests
import random
import json
import time

# Define the GraphQL endpoint
gnomad_api_url = "https://gnomad.broadinstitute.org/api"

# GraphQL query template
REGION_QUERY_TEMPLATE = """
query {
  region(chrom: "%s", start: %d, stop: %d, reference_genome: GRCh38) {
    chrom
    variants(dataset: gnomad_r4) {
      variant_id
      hgvsc
      genome {
        ac
        an
        populations {
          id
          ac
          an
        }
      }
      exome {
        ac
        an
        populations {
          id
          ac
          an
        }
      }
    }
  }
}
"""

# Function to fetch data from gnomAD API
def fetch_data_from_gnomad(chrom, start, stop, retries=3, delay=2):
    chrom = chrom.replace('chr', '')  # Remove 'chr' prefix if present
    query = REGION_QUERY_TEMPLATE % (chrom, start, stop)
    
    for attempt in range(retries):
        response = requests.post(gnomad_api_url, json={"query": query})
        
        if response.status_code == 200:
            data = response.json()
            if data.get("data", {}).get("region", {}).get("variants"):
                return data
        else:
            print(f"Error fetching {chrom}:{start}-{stop}, status: {response.status_code}. Retrying...")
            time.sleep(delay)
    
    return None

# Read the CSV file
csv_file = 'data/SNP_RNA_GC.csv'
df = pd.read_csv(csv_file)

# Select specific rows
first_20_rows = df.head(20)
random_rows = df.iloc[123:].sample(n=100)
selected_rows = pd.concat([first_20_rows, random_rows])

# Prepare results list
results = []
total_rows = len(selected_rows)

# Define population IDs
population_ids = ["remaining", "amr", "fin", "ami", "eas", "mid", "sas", "asj", "afr", "nfe"]

# Loop through selected rows and fetch data
for idx, (_, row) in enumerate(selected_rows.iterrows(), start=1):
    chrom = row['Chromosome']
    start = int(row['Start'])
    stop = int(row['End'])
    gene_name = row['GeneName']

    print(f"Progress: {idx}/{total_rows} - Processing {chrom}:{start}-{stop}")

    # Fetch data from gnomAD
    result = fetch_data_from_gnomad(chrom, start, stop)

    if result:
        variants = result.get('data', {}).get('region', {}).get('variants', [])
        if not variants:
            print(f"No variants found for {chrom}:{start}-{stop}")
        
        for variant in variants:
            variant_id = variant.get('variant_id')
            hgvsc = variant.get('hgvsc')
            
            genome_data = variant.get('genome')
            exome_data = variant.get('exome')
            data_source = genome_data if genome_data else exome_data

            if not data_source:
                continue
            
            ac = data_source.get('ac', 0)
            an = data_source.get('an', 1)
            allele_freq = ac / an if an else 0

            variant_info = {
                'gene_name': gene_name,
                'chrom': chrom,
                'start': start,
                'stop': stop,
                'variation_id': variant_id,
                'variation_consequence': hgvsc,
                'allele_count': ac,
                'allele_num': an,
                'allele_freq': allele_freq,
            }

            # Add population data
            pop_data = {pop["id"]: pop for pop in data_source.get("populations", [])}
            for pop in population_ids:
                variant_info[f"{pop}_ac"] = pop_data.get(pop, {}).get("ac", 0)
                variant_info[f"{pop}_an"] = pop_data.get(pop, {}).get("an", 1)
                variant_info[f"{pop}_af"] = (variant_info[f"{pop}_ac"] / variant_info[f"{pop}_an"]) if variant_info[f"{pop}_an"] else 0
            
            results.append(variant_info)

# Save results to CSV
output_df = pd.DataFrame(results)
output_csv_file = 'data/gnomad_gene_rand_freq.csv'
output_df.to_csv(output_csv_file, index=False)

print(f"Data has been saved to {output_csv_file}")

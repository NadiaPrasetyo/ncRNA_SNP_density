import csv
import requests

gnomad_url = "https://gnomad.broadinstitute.org/api"

gene_query_template = """
query GeneQuery($gene_symbol: String!) {
  gene(gene_symbol: $gene_symbol, reference_genome: GRCh38) {
    chrom
    start
    stop
    hgnc_id
    symbol
    variants(dataset: gnomad_r4) {
      variant_id
      hgvsc
      {DATA_TYPE} {
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

def fetch_gene_data(gene_symbol, data_type="genome"):
    print(f"Fetching {data_type} data for gene: {gene_symbol}")
    query = gene_query_template.replace("{DATA_TYPE}", data_type)
    variables = {"gene_symbol": gene_symbol}
    response = requests.post(gnomad_url, json={"query": query, "variables": variables})
    if response.status_code == 200:
        data = response.json()
        if data.get("data", {}).get("gene") and data["data"]["gene"].get("variants"):
            print(f"Successfully fetched {data_type} data for {gene_symbol}")
            return data
        elif data_type == "genome":
            print(f"No genome data for {gene_symbol}, trying exome...")
            return fetch_gene_data(gene_symbol, data_type="exome")
    print(f"Error fetching data for {gene_symbol}: {response.status_code}")
    return None

def write_results_to_csv(results, output_filepath):
    print(f"Writing results to {output_filepath}")
    population_ids = ["remaining", "amr", "fin", "ami", "eas", "mid", "sas", "asj", "afr", "nfe"]
    fieldnames = [
        "gene_name", "variation_id", "variation_consequence", "allele_count", "allele_num", "allele_freq"
    ] + [f"{pop}_ac" for pop in population_ids] + [f"{pop}_an" for pop in population_ids] + [f"{pop}_af" for pop in population_ids]
    
    with open(output_filepath, "w", newline='', encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for gene, data in results.items():
            if "data" in data and "gene" in data["data"] and data["data"]["gene"]:
                gene_data = data["data"]["gene"]
                for variant in gene_data.get("variants", []):
                    variant_data = variant.get("genome", {}) or variant.get("exome", {})
                    if not variant_data:
                        print(f"Warning: No genome or exome data for {gene} - {variant.get('variant_id', 'Unknown variant')}")
                        continue
                    
                    row = {
                        "gene_name": gene,
                        "variation_id": variant.get("variant_id", ""),
                        "variation_consequence": variant.get("hgvsc", ""),
                        "allele_count": variant_data.get("ac", 0),
                        "allele_num": variant_data.get("an", 1),
                        "allele_freq": variant_data.get("ac", 0) / variant_data.get("an", 1) if variant_data.get("an", 0) else 0
                    }
                    
                    pop_data = {pop["id"]: pop for pop in variant_data.get("populations", [])}
                    for pop in population_ids:
                        row[f"{pop}_ac"] = pop_data.get(pop, {}).get("ac", 0)
                        row[f"{pop}_an"] = pop_data.get(pop, {}).get("an", 1)
                        row[f"{pop}_af"] = pop_data.get(pop, {}).get("ac", 0) / pop_data.get(pop, {}).get("an", 1) if pop_data.get(pop, {}).get("an", 0) else 0
                    
                    writer.writerow(row)
    print("CSV writing complete.")

def main():
    input_filepath = "data/SNP-densities-and-RNA.csv"
    output_filepath = "data/gnomad_gene_data.csv"
    results = {}
    
    with open(input_filepath, newline='', encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        gene_names = [row["GeneName"].strip() for row in reader]
    
    print(f"Processing {len(gene_names)} genes...")
    missing_genes = []
    
    for idx, gene in enumerate(gene_names, start=1):
        print(f"Processing {idx}/{len(gene_names)}: {gene}")
        data = fetch_gene_data(gene)
        if data and data.get("data", {}).get("gene"):
            results[gene] = data
        else:
            print(f"Warning: No data found for gene {gene}")
            missing_genes.append(gene)
    
    write_results_to_csv(results, output_filepath)
    
    if missing_genes:
        print("\nThe following genes had no data:")
        for gene in missing_genes:
            print(gene)

if __name__ == "__main__":
    main()

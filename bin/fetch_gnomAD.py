import csv
import requests
import time

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
REGION_QUERY_TEMPLATE = """
query RegionQuery($chrom: String!, $start: Int!, $stop: Int!) {
  region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
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

def format_chromosome(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom

def fetch_region_data(chrom, start, stop, data_type="genome", retries=3, delay=2):
    chrom = format_chromosome(chrom)
    query = REGION_QUERY_TEMPLATE.replace("{DATA_TYPE}", data_type)
    variables = {"chrom": chrom, "start": start, "stop": stop}
    
    for attempt in range(retries):
        response = requests.post(GNOMAD_URL, json={"query": query, "variables": variables})
        if response.status_code == 200:
            data = response.json()
            if data.get("data", {}).get("region", {}).get("variants"):
                return data
            elif data_type == "genome":
                return fetch_region_data(chrom, start, stop, data_type="exome", retries=retries, delay=delay)
        else:
            print(f"Error fetching {chrom}:{start}-{stop}, status: {response.status_code}. Retrying...")
            time.sleep(delay)
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
        
        for (chrom, start, end, gene_name), data in results.items():
            if "data" in data and "region" in data["data"] and data["data"]["region"]:
                for variant in data["data"]["region"].get("variants", []):
                    variant_data = variant.get("genome", {}) or variant.get("exome", {})
                    if not variant_data:
                        continue
                    
                    row = {
                        "gene_name": gene_name,
                        "variation_id": variant.get("variant_id", ""),
                        "variation_consequence": variant.get("hgvsc", ""),
                        "allele_count": variant_data.get("ac", 0),
                        "allele_num": variant_data.get("an", 1),
                        "allele_freq": variant_data.get("ac", 0) / variant_data.get("an", 1) if variant_data.get("an", 1) else 0
                    }
                    
                    pop_data = {pop["id"]: pop for pop in variant_data.get("populations", [])}
                    for pop in population_ids:
                        row[f"{pop}_ac"] = pop_data.get(pop, {}).get("ac", 0)
                        row[f"{pop}_an"] = pop_data.get(pop, {}).get("an", 1)
                        row[f"{pop}_af"] = row[f"{pop}_ac"] / row[f"{pop}_an"] if row[f"{pop}_an"] else 0
                    
                    writer.writerow(row)
    print("CSV writing complete.")

def main():
    input_filepath = "data/SNP-densities-and-RNA.csv"
    output_filepath = "data/gnomad_region_data.csv"
    results = {}
    
    with open(input_filepath, newline='', encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        regions = [(row["Chromosome"].strip(), int(row["Start"]), int(row["End"]), row["GeneName"].strip()) for row in reader]
    
    missing_regions = []
    for idx, (chrom, start, end, gene_name) in enumerate(regions, start=1):
        print(f"Processing {idx}/{len(regions)}: {chrom}:{start}-{end} ({gene_name})")
        data = fetch_region_data(chrom, start, end)
        if data and data.get("data", {}).get("region"):
            results[(chrom, start, end, gene_name)] = data
        else:
            missing_regions.append((chrom, start, end, gene_name))
    
    write_results_to_csv(results, output_filepath)
    
    if missing_regions:
        print("\nThe following regions had no data:")
        for region in missing_regions:
            print(f"{region[0]}:{region[1]}-{region[2]} ({region[3]})")

if __name__ == "__main__":
    main()

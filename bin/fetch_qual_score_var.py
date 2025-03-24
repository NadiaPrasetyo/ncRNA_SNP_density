import csv
import requests
import time

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"

def format_chromosome(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom

def fetch_variants(chrom, start, end, retries=3, delay=2):
    chrom = format_chromosome(chrom)
    print(f"Fetching variants for {chrom}:{start}-{end}")
    
    query = """
    query GetVariants($chrom: String!, $start: Int!, $end: Int!) {
      region(chrom: $chrom, start: $start, stop: $end, reference_genome: GRCh38) {
        variants(dataset: gnomad_r4) {
          variant_id
        }
      }
    }
    """
    variables = {"chrom": chrom, "start": start, "end": end}

    for attempt in range(retries):
        response = requests.post(GNOMAD_URL, json={'query': query, 'variables': variables})
        if response.status_code == 200:
            try:
                data = response.json()
                if "data" in data and data["data"].get("region", {}).get("variants"):
                    return data["data"]["region"]["variants"]
            except Exception as e:
                print(f"Error parsing JSON for {chrom}:{start}-{end}: {e}")
        else:
            print(f"Request failed ({response.status_code}: {response.text}), retrying in {delay}s...")
            time.sleep(delay)

    print(f"Failed to fetch variants for {chrom}:{start}-{end} after {retries} attempts.")
    return []

def fetch_quality_metrics(variant_id):
    print(f"Fetching quality metrics for variant: {variant_id}")
    query = """
    query GetVariantMetrics($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
        coverage {
          exome { mean }
          genome { mean }
        }
        exome {
          af
          faf95 { popmax popmax_population }
          quality_metrics { site_quality_metrics { metric value } }
        }
        genome {
          af
          faf95 { popmax popmax_population }
          quality_metrics { site_quality_metrics { metric value } }
        }
      }
    }
    """
    variables = {"variantId": variant_id}
    response = requests.post(GNOMAD_URL, json={'query': query, 'variables': variables})

    if response.status_code != 200:
        print(f"Error {response.status_code}: {response.text}")
        return {}

    data = response.json().get("data", {}).get("variant", {})

    if not data:
        print(f"No data found for variant: {variant_id}")
        return {}

    exome_cov = data.get("coverage", {}).get("exome", {}).get("mean")
    genome_cov = data.get("coverage", {}).get("genome", {}).get("mean")

    if exome_cov is None and genome_cov is None:
        print(f"Skipping variant {variant_id} due to lack of coverage data.")
        return {}
    
    genome_cov = genome_cov or 0
    exome_cov = exome_cov or 0
    selected_source = "genome" if genome_cov > exome_cov else "exome"
    selected_data = data.get("genome") if genome_cov > exome_cov else data.get("exome")

    if not selected_data:
        print(f"No quality metrics found for variant: {variant_id}")
        return {}

    quality_metrics = {m["metric"]: m["value"] for m in selected_data.get("quality_metrics", {}).get("site_quality_metrics", [])}

    return {
        "Mean_Coverage": max(genome_cov, exome_cov),
        "seq_data": selected_source,
        "Allele-frequency": selected_data.get("af"),
        "popmax": selected_data.get("faf95", {}).get("popmax_population"),
        "popmax_af": selected_data.get("faf95", {}).get("popmax"),
        **quality_metrics
    }

def main():
    input_file = "data/SNP-densities-and-RNA.csv"
    output_file = "data/variant_qual_metrics.csv"
    
    with open(input_file, newline='') as csvfile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(csvfile)
        fieldnames = [
            "GeneName", "Variation_ID", "Mean_Coverage", "seq_data", "Allele-frequency",
            "popmax", "popmax_af", "SiteQuality", "AS_MQ", "AS_FS", "AS_MQRankSum",
            "AS_pab_max", "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for row in reader:
            chrom, start, end, gene_name = row["Chromosome"], int(row["Start"]), int(row["End"]), row["GeneName"]
            print(f"Processing gene: {gene_name} ({chrom}:{start}-{end})")
            variants = fetch_variants(chrom, start, end)

            for variant in variants:
                variant_id = variant["variant_id"]
                metrics = fetch_quality_metrics(variant_id)

                if not metrics:
                    continue

                writer.writerow({
                    "GeneName": gene_name,
                    "Variation_ID": variant_id,
                    "Mean_Coverage": metrics.get("Mean_Coverage", "N/A"),
                    "seq_data": metrics.get("seq_data", "N/A"),
                    "Allele-frequency": metrics.get("Allele-frequency", "N/A"),
                    "popmax": metrics.get("popmax", "N/A"),
                    "popmax_af": metrics.get("popmax_af", "N/A"),
                    "SiteQuality": metrics.get("SiteQuality", "N/A"),
                    "AS_MQ": metrics.get("AS_MQ", "N/A"),
                    "AS_FS": metrics.get("AS_FS", "N/A"),
                    "AS_MQRankSum": metrics.get("AS_MQRankSum", "N/A"),
                    "AS_pab_max": metrics.get("AS_pab_max", "N/A"),
                    "AS_ReadPosRankSum": metrics.get("AS_ReadPosRankSum", "N/A"),
                    "AS_SOR": metrics.get("AS_SOR", "N/A"),
                    "AS_VarDP": metrics.get("AS_VarDP", "N/A")
                })

                time.sleep(1)  # Prevents hitting API rate limits
    
    print(f"Processing complete. Results saved in {output_file}")

if __name__ == "__main__":
    main()

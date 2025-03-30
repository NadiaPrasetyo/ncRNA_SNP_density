import csv
import requests
import time
import random

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"


def format_chromosome(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom


def is_snp(variant_id):
    parts = variant_id.split("-")
    return len(parts) == 4 and len(parts[2]) == 1 and len(parts[3]) == 1


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
                return data.get("data", {}).get("region", {}).get("variants", [])
            except Exception as e:
                print(f"Error parsing JSON: {e}")
        else:
            print(f"Request failed ({response.status_code}), retrying...")
            time.sleep(delay)
    return []


def fetch_quality_metrics(variant_id):
    if not is_snp(variant_id):
        print(f"Skipping non-SNP variant: {variant_id}")
        return {}

    print(f"Fetching quality metrics for variant: {variant_id}")
    query = """
    query GetVariantMetrics($variantId: String!) {
      variant(variantId: $variantId, dataset: gnomad_r4) {
        coverage { exome { mean } genome { mean } }
        exome { af faf95 { popmax popmax_population } quality_metrics { site_quality_metrics { metric value } } }
        genome { af faf95 { popmax popmax_population } quality_metrics { site_quality_metrics { metric value } } }
      }
    }
    """
    response = requests.post(GNOMAD_URL, json={'query': query, 'variables': {'variantId': variant_id}})

    if response.status_code != 200:
        print(f"Error {response.status_code}")
        return {}

    data = response.json().get("data", {}).get("variant", {})
    exome_cov = data.get("coverage", {}).get("exome", {}).get("mean", 0)
    genome_cov = data.get("coverage", {}).get("genome", {}).get("mean", 0)
    selected_source = "genome" if genome_cov > exome_cov else "exome"
    selected_data = data.get(selected_source, {})

    quality_metrics = {m["metric"]: m["value"] for m in selected_data.get("quality_metrics", {}).get("site_quality_metrics", [])}

    return {
        "Mean_Coverage": max(genome_cov, exome_cov),
        "seq_data": selected_source,
        "Allele-frequency": selected_data.get("af"),
        "popmax": selected_data.get("faf95", {}).get("popmax_population"),
        "popmax_af": selected_data.get("faf95", {}).get("popmax"),
        **quality_metrics
    }

# New functionality: Select 120 random genes from SNP_RNA_GC.csv after line 124
def select_random_genes(input_file):
    with open(input_file, newline='') as csvfile:
        reader = list(csv.DictReader(csvfile))
        filtered_genes = reader[124:]
        return random.sample(filtered_genes, min(120, len(filtered_genes)))


def main():
    input_file = "data/SNP_RNA_GC.csv"
    output_file = "data/variant_qual_metrics_temp.csv"

    selected_genes = select_random_genes(input_file)

    with open(output_file, 'w', newline='') as outfile:
        fieldnames = [
            "GeneName", "Variation_ID", "Mean_Coverage", "seq_data", "Allele-frequency",
            "popmax", "popmax_af", "SiteQuality", "AS_MQ", "AS_FS", "AS_MQRankSum",
            "AS_pab_max", "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP"
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in selected_genes:
            gene_name = row["GeneName"]
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
            variants = fetch_variants(chrom, start, end)

            for variant in variants:
                metrics = fetch_quality_metrics(variant["variant_id"])
                if metrics:
                    writer.writerow({
                        "GeneName": gene_name,
                        "Variation_ID": variant["variant_id"],
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
                time.sleep(1)

    print(f"Processing complete. Results saved in {output_file}")

if __name__ == "__main__":
    main()
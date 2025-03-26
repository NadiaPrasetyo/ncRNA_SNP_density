import csv
import requests
import time

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"

def format_chromosome(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom

def is_snp(variant_id):
    parts = variant_id.split("-")
    if len(parts) == 4 and len(parts[2]) == 1 and len(parts[3]) == 1:
        return True
    return False

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
    if not is_snp(variant_id):
        print(f"Skipping non-SNP variant: {variant_id}")
        return {}
    
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
    skipped_genes = {"GeneX", "GeneY"}  # Add genes to be skipped here
    
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
            gene_name = row["GeneName"]
            if gene_name in skipped_genes:
                print(f"Skipping gene: {gene_name}")
                continue
            
            chrom, start, end = row["Chromosome"], int(row["Start"]), int(row["End"])
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

    
# "RNU5E-1","RNU5D-1","SNORD118","SNAR-B2","TRG-CCC5-1","MIR4538 ","SNAR-C4","RNU5F-1","TRK-TTT3-2","TRE-TTC3-1","TRG1","TRA-CGC5-1-001","tRNA-Gly-CCC-4-1","tRNA-Glu-CTC-1-2","tRNA-Ile-AAT-12-1","tRNA-Gly-GCC-2-4","tRNA-Pro-AGG-2-4","TRG-CCC6-1","TRA-TGC4-1","TRP-AGG1-1","TRK-CTT3-1","tRNA-Glu-TTC-4-1","tRNA-Thr-CGT-5-1","TRG-TCC1-1","SNORD13","tRNA-Val-CAC-3-1","tRNA-Thr-TGT-4-1","TRV-TAC1-1","SNAR-A1","LINC00459","tRNA-Pro-TGG-1-1","FAM30A","tRNA-Gln-CTG-2-1","tRNA-Arg-CCT-3-1","TRV","TRF-GAA1-4","tRNA-iMet-CAT-2-1","SNAR-C3","TRA-AGC11-1","TRL-TAG2-1","tRNA-Met-CAT-6-1","TRM-CAT3-2","tRNA-Leu-CAG-2-2","TRA-AGC1-1","tRNA-Gly-TCC-3-1","tRNA-Lys-CTT-5-1","tRNA-Ile-AAT-2-1","TRA-CGC4-1","TRK-TTT4-1","SNAR-G2","tRNA-Asn-GTT-1-1","TRE-TTC2-1","tRNA-Leu-TAG-3-1","Val-tRNA","TRA-AGC20-1","tRNA-Pro-TGG-3-3","tRNA-Arg-ACG-1-2","TRM-CAT2-1","tRNA-Thr-CGT-6-1","TRA-AGC5-1","tRNA-Val-CAC-2-1","tRNA-Lys-CTT-1-2","TRI-AAT4-1","tRNA-Lys-TTT-1-1","TRA-AGC4-1","tRNA-Lys-CTT-2-1","TRI-AAT5-1","tRNA-Thr-TGT-1-1","tRNA-Asp-GTC-1-1","tRNA-Lys-TTT-6-1","TRF-GAA2-1","TRP-TGG2-1","TRT-TGT5-1","TRR-TCG1-1","MIR3689A","tRNA-Ala-CGC-3-1","TRE-TTC5-1","TRV-AAC5-1","TRR-ACG2-3","TRT-AGT1-2","TRL-TAA4-1","tRNA-Ala-TGC-3-1","SYT15-AS1","TRR-TCG3-1","TRS-AGA2-6","tRNA-Ser-TGA-2-1","lnc-SLCO4A1-8","tRNA-Ile-AAT-1-1","tRNA-Gln-TTG-3-1","TRQ-CTG6-1","TRK-TTT7-1","TRT-TGT2-1","TRV-CAC1-5","tRNA-Ala-CGC-1-1","TRX-CAT1-5","tRNA-Val-CAC-1-3","SNAR-G1","tRNA-Arg-TCG-4-1","TRC-GCA11-1","TRD-GTC2-8","TRA-AGC15-1","tRNA-Ser-GCT-2-1","tRNA-Ser","RNU6-9","TRC-GCA4-1","TRL-AAG2-4","TRN-GTT2-4","TRA-CGC2-1","RCCD1-AS1","TRG-CCC2-1","TRQ-TTG2-1","tRNA-Val-AAC-3-1","tRNA-Ala-TGC-2-1","tRNA-Asn-GTT-4-1","tRNA-Ser-CGA-2-1","tRNA-Asn-GTT-3-2","LINC01671","SNAR-F","BLACE","tRNA-Lys-TTT-5-1","tRNA-Arg-CCT-2-1","TRS-AGA2-6"
import requests
import csv
import time

def query_gnomad(chrom, start, end, gene_name):
    url = "https://gnomad.broadinstitute.org/api"
    query = """
    query VariantsInRegion($chrom: String!, $start: Int!, $end: Int!, $reference_genome: ReferenceGenomeId!, $dataset: DatasetId!) {
      region(chrom: $chrom, start: $start, stop: $end, reference_genome: $reference_genome) {
        variants(dataset: $dataset) {
          variant_id
          genome {
            af
            populations {
              id
              ac
              an
              af
            }
          }
          exome {
            af
            populations {
              id
              ac
              an
              af
            }
          }
        }
      }
    }
    """
    variables = {
        "chrom": chrom,
        "start": start,
        "end": end,
        "reference_genome": "GRCh38",  # Use the correct genome version
        "dataset": "gnomad_r4"  # Correctly pass the dataset inside the "variants" field
    }
    
    response = requests.post(url, json={"query": query, "variables": variables})
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error querying gnomAD: {response.text}")
        return None

def main():
    input_file = "data/SNP-densities-and-RNA.csv"
    output_file = "gnomad_variant_frequencies.csv"
    
    with open(input_file, newline='') as csvfile, open(output_file, mode='w', newline='') as outfile:
        reader = csv.DictReader(csvfile)
        fieldnames = ["Chromosome", "Start", "End", "GeneName", "VariantID", 
                      "Genome_AF", "Exome_AF", "Population", "AC", "AN", "AF"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for row in reader:
            chrom = row["Chromosome"].replace("chr", "")  # Ensure chromosome formatting
            start = int(row["Start"])
            end = int(row["End"])
            gene_name = row["GeneName"]
            
            print(f"Querying variants for {gene_name} (chr{chrom}:{start}-{end})")
            result = query_gnomad(chrom, start, end, gene_name)
            
            if result and "data" in result and result["data"].get("region"):
                for variant in result["data"]["region"]["variants"]:
                    variant_id = variant["variant_id"]
                    genome_af = variant["genome"]["af"] if variant.get("genome") else "NA"
                    exome_af = variant["exome"]["af"] if variant.get("exome") else "NA"
                    
                    # Extract population frequency data
                    populations = []
                    if "genome" in variant and variant["genome"]:
                        populations.extend(variant["genome"].get("populations", []))
                    if "exome" in variant and variant["exome"]:
                        populations.extend(variant["exome"].get("populations", []))
                    
                    for pop in populations:
                        writer.writerow({
                            "Chromosome": chrom,
                            "Start": start,
                            "End": end,
                            "GeneName": gene_name,
                            "VariantID": variant_id,
                            "Genome_AF": genome_af,
                            "Exome_AF": exome_af,
                            "Population": pop["id"],
                            "AC": pop["ac"],
                            "AN": pop["an"],
                            "AF": pop["af"]
                        })
            else:
                print(f"No variants found for {gene_name}.")
            
            time.sleep(1)  # Avoid excessive API requests

if __name__ == "__main__":
    main()

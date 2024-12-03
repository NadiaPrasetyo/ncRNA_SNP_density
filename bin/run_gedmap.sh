#!/bin/bash

# Static file paths
FASTA_FILE="data/Chromosomes_FASTA/hg38_chromosomes.fa"
VCF_FILE="data/filtered_variants.vcf"
OUTPUT_DIR="data/GEDS"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Read the CSV file and process each line (skip the header row)
while IFS=, read -r chromosome start end gene_id gene_name snp_density rna_type; do
    # Skip header row
    if [[ "$gene_name" == "GeneName" ]]; then
        continue
    fi
    
    # Define the output prefix for this gene
    output_prefix="$OUTPUT_DIR/$gene_name"

    echo "Running gedmap parse for gene: $gene_name"
    ../gedmap/gedmap parse "$FASTA_FILE" "$VCF_FILE" -o "$output_prefix"
done < data/SNP-densities-and-RNA.csv

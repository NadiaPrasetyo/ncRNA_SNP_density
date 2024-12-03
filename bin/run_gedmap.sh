#!/bin/bash

# Directory paths
FASTA_DIR="data/FASTA"
VCF_DIR="data/VCF"
OUTPUT_DIR="data/GEDS"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Read the CSV file and process each line (skip the header row)
while IFS=, read -r chromosome start end gene_id gene_name snp_density rna_type; do
    # Skip header row
    if [[ "$gene_name" == "GeneName" ]]; then
        continue
    fi

    #skip the two genes that are HUGE
    if [[ "$gene_name" == "TRA-AGC5-1" || "$gene_name" == "TRE-TTC5-1" ]]; then
        continue
    fi

    # Define the FASTA and VCF file paths based on GeneName
    fasta_file="$FASTA_DIR/$gene_name.fa"
    vcf_file="$VCF_DIR/${gene_name}_variants.vcf"

    # Check if both FASTA and VCF files exist
    if [[ -f "$fasta_file" && -f "$vcf_file" ]]; then
        # Run gedmap parse for the gene
        output_prefix="$OUTPUT_DIR/$gene_name"
        echo "Running gedmap parse for gene: $gene_name"
        ../gedmap/gedmap parse "$fasta_file" "$vcf_file" -o "$output_prefix"
    else
        echo "fasta file: $fasta_file"
        echo "vcf file: $vcf_file"
        echo "Skipping $gene_name: FASTA or VCF file not found."
    fi
done < data/SNP-densities-and-RNA.csv

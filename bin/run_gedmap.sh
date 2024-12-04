#!/bin/bash

# Directory paths
FASTA_DIR="data/Chromosomes_FASTA"
VCF_DIR="data/VCF"
OUTPUT_DIR="data/GEDS"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through the chromosomes (chr1 to chr22)
for i in {1..22}; do
    chromosome="chr$i"
    
    # Define the FASTA and VCF file paths based on chromosome
    fasta_file="$FASTA_DIR/${chromosome}.fa"
    vcf_file="$VCF_DIR/${chromosome}.vcf"

    # Check if both FASTA and VCF files exist
    if [[ -f "$fasta_file" && -f "$vcf_file" ]]; then
        # Run gedmap parse for the chromosome
        output_prefix="$OUTPUT_DIR/$chromosome"
        echo "Running gedmap parse for chromosome: $chromosome"
        ../gedmap/gedmap parse "$fasta_file" "$vcf_file" -o "$output_prefix"
    else
        echo "fasta file: $fasta_file"
        echo "vcf file: $vcf_file"
        echo "Skipping $chromosome: FASTA or VCF file not found."
    fi
done

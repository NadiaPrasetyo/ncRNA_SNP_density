#!/bin/bash

# Script to process FASTA and VCF files for chromosomes 1 to 22 using gedmap.
# It checks for the existence of required files and outputs results to a specified directory.

# Exit immediately if a command exits with a non-zero status.
set -e

# Enable debugging for more detailed error reporting (optional).
set -x

# Directory paths
FASTA_DIR="data/Chromosomes_FASTA"  # Directory containing FASTA files
VCF_DIR="data/VCF"                  # Directory containing VCF files
OUTPUT_DIR="data/GEDS"              # Directory to store output

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through the chromosomes (chr1 to chr22)
for i in {1..22}; do
    chromosome="chr$i"

    # Define the FASTA and VCF file paths based on chromosome
    fasta_file="$FASTA_DIR/${chromosome}.fa"
    vcf_file="$VCF_DIR/${chromosome}.vcf"

    # Check if both FASTA and VCF files exist
    if [[ -f "$fasta_file" && -f "$vcf_file" ]]; then
        output_prefix="$OUTPUT_DIR/$chromosome"  # Prefix for output files
        echo "Running gedmap parse for chromosome: $chromosome"

        # Execute gedmap parse, catching errors
        if ! ../gedmap/gedmap parse "$fasta_file" "$vcf_file" -o "$output_prefix"; then
            echo "Error: gedmap parse failed for chromosome $chromosome" >&2
            continue
        fi
    else
        # Report missing files
        if [[ ! -f "$fasta_file" ]]; then
            echo "Missing FASTA file: $fasta_file" >&2
        fi
        if [[ ! -f "$vcf_file" ]]; then
            echo "Missing VCF file: $vcf_file" >&2
        fi
        echo "Skipping $chromosome: Required files not found."
    fi
done

# End of script

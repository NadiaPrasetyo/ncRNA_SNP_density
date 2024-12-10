#!/bin/bash

# Script to process a single FASTA file (hg38.fa) with all VCF files in a directory using gedmap.
# It checks for the existence of the required files and outputs results to a specified directory.

# Exit immediately if a command exits with a non-zero status.
set -e

# Enable debugging for more detailed error reporting (optional).
set -x

# File and directory paths
FASTA_FILE="data/hg38.fa"  # Path to the single FASTA file
VCF_DIR="data/VCF/processed"         # Directory containing VCF files
OUTPUT_DIR="data/GEDS/genomes"  # Directory to store output

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the FASTA file exists
if [[ ! -f "$FASTA_FILE" ]]; then
    echo "Error: FASTA file $FASTA_FILE not found." >&2
    exit 1
fi

# Loop through all VCF files in the directory
for vcf_file in "$VCF_DIR"/*.vcf; do
    # Check if the VCF file exists
    if [[ -f "$vcf_file" ]]; then
        # Extract the base name of the VCF file for output naming
        base_name=$(basename "$vcf_file" .vcf)
        output_prefix="$OUTPUT_DIR/$base_name"  # Prefix for output files

        # Check if output files already exist
        if [[ -f "${output_prefix}.geds" ]]; then
            echo "Output file ${output_prefix}.geds already exists. Skipping $vcf_file."
            continue
        fi

        echo "Running gedmap parse for VCF file: $vcf_file"

        # Execute gedmap parse, catching errors
        if ! ../gedmap/gedmap parse "$FASTA_FILE" "$vcf_file" -o "$output_prefix"; then
            echo "Error: gedmap parse failed for VCF file $vcf_file" >&2
            continue
        fi
    else
        echo "No VCF files found in $VCF_DIR. Skipping."
        break
    fi
done

# End of script

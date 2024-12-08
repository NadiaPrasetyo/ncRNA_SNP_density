#!/bin/bash

# Script to process GEDS, create minimizer indices, and align FASTQ files for chromosomes 1 to 22.
# Additionally, processes genes from a CSV file for SNP densities and RNA alignment.

# Exit immediately if a command exits with a non-zero status.
set -e

# Enable debugging for detailed command execution logs (optional).
# set -x

# Directory paths
OUTPUT_DIR="data/GEDS"               # Directory for GEDS files
ALIGN_OUTPUT_DIR="data/Alignments"   # Directory for alignment outputs
FASTQ_DIR="data/fastq"               # Directory containing FASTQ files
CSV_FILE="data/SNP-densities-and-RNA.csv"  # CSV file containing gene information

# Chromosomes list
CHROMOSOMES=()
for i in {1..22}; do
    CHROMOSOMES+=("chr$i")
done

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$FASTQ_DIR"

# Step 1: Process each chromosome
for CHROM in "${CHROMOSOMES[@]}"; do
    # Define file paths for GEDS and static files
    GEDS_FILE="$OUTPUT_DIR/$CHROM.geds"
    GEDMAP_POS2FA_FILE="$OUTPUT_DIR/$CHROM.geds.2fa"
    GEDMAP_INDEX_FILE="$OUTPUT_DIR/$CHROM.geds.min"

    # Check for required files
    if [[ ! -f "$GEDS_FILE" || ! -f "$GEDMAP_POS2FA_FILE" ]]; then
        echo "Error: Static GEDS or 2FA file not found for $CHROM. Ensure $GEDS_FILE and $GEDMAP_POS2FA_FILE are present." >&2
        exit 1
    fi

    # Check or create minimizer index file
    if [[ ! -f "$GEDMAP_INDEX_FILE" ]]; then
        echo "Running GEDMAP index for $CHROM..."
        if ! ../gedmap/gedmap index "$GEDS_FILE" -2fa "$GEDMAP_POS2FA_FILE" -o "$GEDMAP_INDEX_FILE" -t 8; then
            echo "Error: GEDMAP index failed for $CHROM." >&2
            exit 1
        fi
    else
        echo "Index file already exists for $CHROM: $GEDMAP_INDEX_FILE"
    fi
done

# Step 2: Process genes from the CSV file
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: CSV file not found: $CSV_FILE" >&2
    exit 1
fi

while IFS=, read -r chromosome start end gene_id gene_name snp_density rna_type; do
    # Skip the header row
    if [[ "$gene_name" == "GeneName" ]]; then
        continue
    fi

    echo "Processing gene: $gene_name"

    for chromosome in "${CHROMOSOMES[@]}"; do
        # Define file paths
        GEDS_FILE="$OUTPUT_DIR/$chromosome.geds"
        GEDMAP_INDEX_FILE="$OUTPUT_DIR/$chromosome.geds.min"
        fastq_file="$FASTQ_DIR/$gene_name.fq"
        sam_file="$ALIGN_OUTPUT_DIR/$gene_name-$chromosome.sam"

        echo " GEDS file: $GEDS_FILE, GEDMAP index file: $GEDMAP_INDEX_FILE, FASTQ file: $fastq_file, SAM file: $sam_file"

        # Check if FASTQ file exists
        if [[ -f "$fastq_file" ]]; then
            # Align FASTQ to GEDS if the SAM file doesn't already exist
            if [[ ! -f "$sam_file" ]]; then
                echo "  Running GEDMAP align for $gene_name on chromosome $chromosome..."
                if ! ../gedmap/gedmap align "$fastq_file" "$GEDS_FILE" "$GEDMAP_INDEX_FILE" -o "$sam_file"; then
                    echo "  Error: GEDMAP align failed for $gene_name on chromosome $chromosome." >&2
                    continue
                fi
            else
                echo "  SAM file already exists for $gene_name on chromosome $chromosome: $sam_file"
            fi

            echo "  Alignment completed for $gene_name on chromosome $chromosome. Output SAM file: $sam_file"
        else
            echo "  Skipping $gene_name: FASTQ file not found."
            echo "  FASTQ file: $fastq_file"
        fi
    done
done < "$CSV_FILE"

echo "Script completed successfully."

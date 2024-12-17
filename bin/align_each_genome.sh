#!/bin/bash

# Script to process GEDS, create minimizer indices, and align FASTQ files for specified genomes.
# Additionally, processes genes from a CSV file for SNP densities and RNA alignment.

# Exit immediately if a command exits with a non-zero status.
set -e

# Enable debugging for detailed command execution logs (optional).
# set -x

# Directory paths
OUTPUT_DIR="data/GEDS/genomes"       # Directory for genome-specific GEDS files
ALIGN_OUTPUT_DIR="data/Alignments"  # Directory for alignment outputs
FASTQ_DIR="data/fastq"              # Directory containing FASTQ files
CSV_FILE="data/SNP-densities-and-RNA.csv"  # CSV file containing gene information

# Genome list
GENOMES=(
    HG00438 HG00621 HG00673 HG00733 HG00735 HG00741 HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258
    HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109 HG02145 HG02148 HG02257 HG02486
    HG02559 HG02572 HG02622 HG02630 HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492 HG03516
    HG03540 HG03579 NA18906 NA20129 NA21309
)

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ALIGN_OUTPUT_DIR"
mkdir -p "$FASTQ_DIR"

# Step 1: Process each genome
for GENOME in "${GENOMES[@]}"; do
    # Define file paths for GEDS and static files
    GEDS_FILE="$OUTPUT_DIR/processed_$GENOME.geds"
    GEDMAP_INDEX_FILE="$OUTPUT_DIR/$GENOME.geds.min"

    # Check for required files
    if [[ ! -f "$GEDS_FILE" ]]; then
        echo "Error: Static GEDS file not found for $GENOME. Ensure $GEDS_FILE are present." >&2
        exit 1
    fi

    # Check or create minimizer index file
    if [[ ! -f "$GEDMAP_INDEX_FILE" ]]; then
        echo "Running GEDMAP index for $GENOME..."
        if ! ../gedmap/gedmap index "$GEDS_FILE" -o "$GEDMAP_INDEX_FILE" -t 8; then
            echo "Error: GEDMAP index failed for $GENOME." >&2
            exit 1
        fi
    else
        echo "Index file already exists for $GENOME: $GEDMAP_INDEX_FILE"
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

    for GENOME in "${GENOMES[@]}"; do
        # Define file paths
        GEDS_FILE="$OUTPUT_DIR/processed_$GENOME.geds"
        GEDMAP_INDEX_FILE="$OUTPUT_DIR/$GENOME.geds.min"
        fastq_file="$FASTQ_DIR/$gene_name.fq"
        sam_file="$ALIGN_OUTPUT_DIR/$gene_name-$GENOME.sam"

        echo " GEDS file: $GEDS_FILE, GEDMAP index file: $GEDMAP_INDEX_FILE, FASTQ file: $fastq_file, SAM file: $sam_file"

        # Check if FASTQ file exists
        if [[ -f "$fastq_file" ]]; then
            # Align FASTQ to GEDS if the SAM file doesn't already exist
            if [[ ! -f "$sam_file" ]]; then
                echo "  Running GEDMAP align for $gene_name on genome $GENOME..."
                if ! ../gedmap/gedmap align "$fastq_file" "$GEDS_FILE" "$GEDMAP_INDEX_FILE" -o "$sam_file"; then
                    echo "  Error: GEDMAP align failed for $gene_name on genome $GENOME." >&2
                    continue
                fi
            else
                echo "  SAM file already exists for $gene_name on genome $GENOME: $sam_file"
            fi

            echo "  Alignment completed for $gene_name on genome $GENOME. Output SAM file: $sam_file"
        else
            echo "  Skipping $gene_name: FASTQ file not found."
            echo "  FASTQ file: $fastq_file"
        fi
    done
done < "$CSV_FILE"

echo "Script completed successfully."

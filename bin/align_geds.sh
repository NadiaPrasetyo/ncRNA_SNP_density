#!/bin/bash

# Directory paths
OUTPUT_DIR="data/GEDS"
ALIGN_OUTPUT_DIR="data/Alignments"  # Directory for alignment outputs
FASTQ_DIR="data/fastq"  # Directory where dynamic fastq files are stored

# Chromosomes list
CHROMOSOMES=()
for i in {1..22}; do
    CHROMOSOMES+=("chr$i")
done


# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $ALIGN_OUTPUT_DIR
mkdir -p $FASTQ_DIR

# Loop through each chromosome
for CHROM in "${CHROMOSOMES[@]}"; do
    # Define file paths for the current chromosome's GEDS and static files
    GEDS_FILE="$OUTPUT_DIR/$CHROM.geds"
    GEDMAP_POS2FA_FILE="$OUTPUT_DIR/$CHROM.geds.2fa"
    GEDMAP_INDEX_FILE="$OUTPUT_DIR/$CHROM.geds.min"

    # Check if the GEDS and static files exist
    if [[ ! -f "$GEDS_FILE" || ! -f "$GEDMAP_POS2FA_FILE" ]]; then
        echo "Error: Static GEDS or 2FA file not found for $CHROM. Please ensure $GEDS_FILE and $GEDMAP_POS2FA_FILE are present."
        exit 1
    fi

    # Check if the minimizer index file exists for this chromosome, if not, create it
    if [[ ! -f "$GEDMAP_INDEX_FILE" ]]; then
        echo "Running GEDMAP index for $CHROM..."
        ../gedmap/gedmap index "$GEDS_FILE" -2fa "$GEDMAP_POS2FA_FILE" -o "$GEDMAP_INDEX_FILE" -t 8

        if [[ $? -ne 0 ]]; then
            echo "Error: GEDMAP index failed for $CHROM."
            exit 1
        fi
    else
        echo "Index file already exists for $CHROM: $GEDMAP_INDEX_FILE"
    fi
done

# Read the CSV file and process each line (skip the header row)
while IFS=, read -r chromosome start end gene_id gene_name snp_density rna_type; do
    # Skip header row
    if [[ "$gene_name" == "GeneName" ]]; then
        continue
    fi

    # Check if the gene corresponds to one of the chromosomes
    if [[ ! " ${CHROMOSOMES[@]} " =~ " $chromosome " ]]; then
        echo "Skipping $gene_name: Chromosome $chromosome not in list."
        continue
    fi

    # Define file paths for the gene-specific FASTQ and SAM output
    fastq_file="$FASTQ_DIR/$gene_name.fq"
    sam_file="$ALIGN_OUTPUT_DIR/$gene_name-$chromosome.sam"

    echo "Processing gene: $gene_name on chromosome $chromosome"

    # Check if the required FASTQ file exists
    if [[ -f "$fastq_file" ]]; then
        # Align: Align FASTQ to GEDS for this chromosome (only if the SAM file doesn't exist)
        if [[ ! -f "$sam_file" ]]; then
            echo "Running GEDMAP align for $gene_name on chromosome $chromosome..."
            ../gedmap/gedmap align "$fastq_file" "$GEDS_FILE" "$GEDMAP_INDEX_FILE" -o "$sam_file" -rc -mao 1 -2fa "$GEDMAP_POS2FA_FILE" -d 10
            if [[ $? -ne 0 ]]; then
                echo "Error: GEDMAP align failed for $gene_name on chromosome $chromosome."
                continue
            fi
        else
            echo "SAM file already exists for $gene_name on chromosome $chromosome: $sam_file"
        fi

        echo "Alignment completed for $gene_name on chromosome $chromosome. Output SAM file: $sam_file"
    else
        echo "Skipping $gene_name: FASTQ file not found."
        echo "FASTQ file: $fastq_file"
    fi
done < data/SNP-densities-and-RNA.csv

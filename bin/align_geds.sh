#!/bin/bash

# Directory paths
FASTQ_FILE="data/ERR187509.fq"  # Static FASTQ file
OUTPUT_DIR="data/GEDS"
ALIGN_OUTPUT_DIR="data/Alignments"  # Directory for alignment outputs

# Create output directories if they don't exist
mkdir -p $OUTPUT_DIR
mkdir -p $ALIGN_OUTPUT_DIR

# Read the CSV file and process each line (skip the header row)
while IFS=, read -r chromosome start end gene_id gene_name snp_density rna_type; do
    # Skip header row
    if [[ "$gene_name" == "GeneName" ]]; then
        continue
    fi

    # Define file paths based on GeneName
    geds_file="$OUTPUT_DIR/$gene_name.geds"
    adj_file="$geds_file.adj"
    pos2fa_file="$geds_file.2fa"
    index_file="$OUTPUT_DIR/$gene_name.min"
    sam_file="$ALIGN_OUTPUT_DIR/$gene_name.sam"

    echo "Processing gene: $gene_name"

    # Check if the required GEDS file exists
    if [[ -f "$geds_file" ]]; then
        # Index: Generate minimizer index
        if [[ ! -f "$index_file" ]]; then
            echo "Running gedmap index for $gene_name..."
            if [[ -f "$adj_file" ]]; then
                ../gedmap/gedmap index "$geds_file" -2fa "$pos2fa_file" -k 3 -w 3 -o "$index_file" -a "$adj_file"
            else
                ../gedmap/gedmap index "$geds_file" -2fa "$pos2fa_file" -k 3 -w 3 -o "$index_file"
            fi
            if [[ $? -ne 0 ]]; then
                echo "Error: GEDMAP index failed for $gene_name."
                continue
            fi
        else
            echo "Index file already exists: $index_file"
        fi

        # Align: Align FASTQ to GEDS
        if [[ ! -f "$sam_file" ]]; then
            echo "Running gedmap align for $gene_name..."
            ../gedmap/gedmap align "$FASTQ_FILE" "$geds_file" "$index_file" -o "$sam_file" -rc -mao 1 -2fa "$pos2fa_file" -d 10
            if [[ $? -ne 0 ]]; then
                echo "Error: GEDMAP align failed for $gene_name."
                continue
            fi
        else
            echo "SAM file already exists: $sam_file"
        fi

        echo "Alignment completed for $gene_name. Output SAM file: $sam_file"
    else
        echo "Skipping $gene_name: GEDS file not found."
        echo "GEDS file: $geds_file"
    fi
done < data/SNP-densities-and-RNA.csv

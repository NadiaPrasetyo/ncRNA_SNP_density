#!/bin/bash

# Static file paths
FASTA_FILE="data/Chromosomes_FASTA/hg38_chromosomes.fa"
VCF_FILE="data/filtered_variants.vcf"
OUTPUT_PREFIX="data/GEDS/hg38"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

../gedmap/gedmap parse "$FASTA_FILE" "$VCF_FILE" -o "$output_prefix"

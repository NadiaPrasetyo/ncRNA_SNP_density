#!/bin/bash

# ----------- CONFIG -----------

# Path to reference genome
REFERENCE="GRCh38.primary_assembly.genome.fa"

# Path to CSV with region definitions
CSV="data/SNP-densities-and-RNA.csv"

# Directory containing FASTQ files (can be in subfolders)
FASTQ_ROOT="data/1000genomes/test/data/NA12282/sequence_read"

# Output directory
OUTPUT_DIR="output"
mkdir -p "$OUTPUT_DIR"

# ----------- PREP -----------

# Convert CSV to BED (remove header)
BED="$OUTPUT_DIR/regions.bed"
tail -n +2 "$CSV" | awk -F',' 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$BED"

# ----------- MAIN LOOP -----------

# Find all FASTQ files recursively
find "$FASTQ_ROOT" -type f -name "*.fastq" | while read -r FASTQ; do

    # Extract filename prefix (without path or extension)
    BASENAME=$(basename "$FASTQ" .fastq)
    PREFIX="${BASENAME%%.*}"  # in case of .fastq.gz or paired files

    echo "Processing $FASTQ → $PREFIX"

    # 1. Align to reference genome
    bwa mem -t 4 "$REFERENCE" "$FASTQ" > "$OUTPUT_DIR/$PREFIX.sam"

    # 2. Convert to sorted BAM
    samtools view -S -b "$OUTPUT_DIR/$PREFIX.sam" | samtools sort -o "$OUTPUT_DIR/$PREFIX.sorted.bam"
    samtools index "$OUTPUT_DIR/$PREFIX.sorted.bam"

    # 3. Filter reads that map to regions
    bedtools intersect -abam "$OUTPUT_DIR/$PREFIX.sorted.bam" -b "$BED" > "$OUTPUT_DIR/$PREFIX.filtered.bam"
    samtools index "$OUTPUT_DIR/$PREFIX.filtered.bam"

    # 4. Count how many reads passed the filter
    READ_COUNT=$(samtools view "$OUTPUT_DIR/$PREFIX.filtered.bam" | wc -l)

    if [[ "$READ_COUNT" -gt 0 ]]; then
        echo "✅ $PREFIX has $READ_COUNT reads mapped to regions."
    else
        echo "❌ $PREFIX has NO reads mapped to regions."
    fi

    echo "---"

done

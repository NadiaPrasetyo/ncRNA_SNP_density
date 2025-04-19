#!/bin/bash

# ----------- CONFIG -----------

REFERENCE="GRCh38.primary_assembly.genome.fa"
CSV="data/SNP-densities-and-RNA.csv"
FASTQ_ROOT="data/1000genomes"
OUTPUT_DIR="data/1000genomes/output"
mkdir -p "$OUTPUT_DIR"

# Final combined BED output
COMBINED_BED="$OUTPUT_DIR/all_filtered_reads.bed"
> "$COMBINED_BED"  # clear if exists

# ----------- PREP -----------

BED="$OUTPUT_DIR/regions.bed"
tail -n +2 "$CSV" | awk -F',' 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$BED"

# ----------- MAIN LOOP -----------

find "$FASTQ_ROOT" -type f -name "*.fastq" | while read -r FASTQ; do

    BASENAME=$(basename "$FASTQ" .fastq)
    PREFIX="${BASENAME%%.*}"  # handle .fastq.gz too

    echo "Processing $FASTQ → $PREFIX"

    # 1. Align
    bwa mem -t 4 "$REFERENCE" "$FASTQ" > "$OUTPUT_DIR/$PREFIX.sam"

    # 2. Sort and index BAM
    samtools view -S -b "$OUTPUT_DIR/$PREFIX.sam" | samtools sort -o "$OUTPUT_DIR/$PREFIX.sorted.bam"
    samtools index "$OUTPUT_DIR/$PREFIX.sorted.bam"

    # 3. Filter to regions
    bedtools intersect -abam "$OUTPUT_DIR/$PREFIX.sorted.bam" -b "$BED" > "$OUTPUT_DIR/$PREFIX.filtered.bam"
    samtools index "$OUTPUT_DIR/$PREFIX.filtered.bam"

    # 4. Count reads
    READ_COUNT=$(samtools view "$OUTPUT_DIR/$PREFIX.filtered.bam" | wc -l)

    if [[ "$READ_COUNT" -gt 0 ]]; then
        echo "✅ $PREFIX has $READ_COUNT reads mapped to regions."

        # 5. Convert filtered reads to BED and tag with sample name
        bedtools bamtobed -i "$OUTPUT_DIR/$PREFIX.filtered.bam" | \
        awk -v sample="$PREFIX" 'BEGIN{OFS="\t"} {print $0, sample}' >> "$COMBINED_BED"
    else
        echo "❌ $PREFIX has NO reads mapped to regions."
    fi

    echo "---"

done

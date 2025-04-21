#!/bin/bash

# ----------- CONFIG -----------

REFERENCE="GRCh38.primary_assembly.genome.fa"
CSV="data/SNP-densities-and-RNA.csv"
FASTQ_ROOT="data/1000genomes"
OUTPUT_DIR="data/1000genomes/output"
mkdir -p "$OUTPUT_DIR"

# Final combined outputs
COMBINED_BED="$OUTPUT_DIR/all_filtered_reads.bed"
COMBINED_BAM="$OUTPUT_DIR/combined.filtered.bam"
TMP_BAM_DIR="$OUTPUT_DIR/tmp_bams"
mkdir -p "$TMP_BAM_DIR"

> "$COMBINED_BED"  # Clear if exists

# ----------- PREP -----------

BED="$OUTPUT_DIR/regions.bed"
tail -n +2 "$CSV" | awk -F',' 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$BED"

# ----------- MAIN LOOP -----------

FILTERED_BAMS=()

find "$FASTQ_ROOT" -type f -name "*.fastq" | while read -r FASTQ; do

    BASENAME=$(basename "$FASTQ" .fastq)
    PREFIX="${BASENAME%%.*}"

    echo "Processing $FASTQ → $PREFIX"

    SORTED_BAM="$TMP_BAM_DIR/$PREFIX.sorted.bam"
    FILTERED_BAM="$TMP_BAM_DIR/$PREFIX.filtered.bam"

    # 1. Align, convert, sort BAM
    bwa mem -t 2 "$REFERENCE" "$FASTQ" | \
    samtools view -Sb - | \
    samtools sort -o "$SORTED_BAM"

    if [ ! -s "$SORTED_BAM" ]; then
        echo "❌ Alignment failed for $PREFIX, skipping..."
        echo "---"
        continue
    fi

    samtools index "$SORTED_BAM"

    # 2. Filter to regions
    bedtools intersect -abam "$SORTED_BAM" -b "$BED" > "$FILTERED_BAM"
    samtools index "$FILTERED_BAM"

    READ_COUNT=$(samtools view "$FILTERED_BAM" | wc -l)

    if [[ "$READ_COUNT" -gt 0 ]]; then
        echo "✅ $PREFIX has $READ_COUNT reads mapped to regions."
        FILTERED_BAMS+=("$FILTERED_BAM")
        bedtools bamtobed -i "$FILTERED_BAM" | \
        awk -v sample="$PREFIX" 'BEGIN{OFS="\t"} {print $0, sample}' >> "$COMBINED_BED"
    else
        echo "❌ $PREFIX has NO reads mapped to regions."
    fi

    echo "---"

done

# ----------- MERGE FILTERED BAMs -----------

if [ "${#FILTERED_BAMS[@]}" -gt 0 ]; then
    echo "Merging ${#FILTERED_BAMS[@]} filtered BAMs..."
    samtools merge -f "$COMBINED_BAM" "${FILTERED_BAMS[@]}"
    samtools index "$COMBINED_BAM"
    echo "✅ Combined BAM created at: $COMBINED_BAM"
else
    echo "❌ No BAMs to merge. Exiting."
fi

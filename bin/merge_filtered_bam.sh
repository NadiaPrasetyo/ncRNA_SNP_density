#!/bin/bash

# ----------- CONFIG -----------
TMP_BAM_DIR="data/1000genomes/output/tmp_bams"
OUTPUT_DIR="data/1000genomes/output"
COMBINED_BAM="$OUTPUT_DIR/combined.filtered.bam"

# ----------- FIND AND MERGE FILTERED BAMs -----------

# Find all filtered BAMs in tmp directory
FILTERED_BAMS=($(find "$TMP_BAM_DIR" -name "*.filtered.bam"))

if [ "${#FILTERED_BAMS[@]}" -eq 0 ]; then
    echo "❌ No filtered BAM files found to merge in $TMP_BAM_DIR"
    exit 1
fi

echo "🔍 Found ${#FILTERED_BAMS[@]} filtered BAM files to merge..."
for bam in "${FILTERED_BAMS[@]}"; do
    echo "➕ $bam"
done

# Merge BAMs
samtools merge -f "$COMBINED_BAM" "${FILTERED_BAMS[@]}"
samtools index "$COMBINED_BAM"

echo "✅ Combined BAM saved as: $COMBINED_BAM"

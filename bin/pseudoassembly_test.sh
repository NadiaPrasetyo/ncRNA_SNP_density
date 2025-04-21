#!/bin/bash

set -euo pipefail

# Set paths
FILTERED_BAM="data/1000genomes/output/combined.filtered.bam"
FILTERED_FASTQ="data/1000genomes/output/combined.filtered.fastq"
ASSEMBLY_DIR="data/1000genomes/output/megahit_filtered_assembly"
REFERENCE="GRCh38.primary_assembly.genome.fa"
CONTIGS_FASTA="$ASSEMBLY_DIR/final.contigs.fa"
CONTIG_ALIGN_SAM="data/1000genomes/output/contigs_vs_reference.sam"

# Step 1: Convert BAM to FASTQ
echo "Converting filtered BAM to FASTQ..."
samtools fastq "$FILTERED_BAM" > "$FILTERED_FASTQ"

# Step 2: Calculate stats
NUM_READS=$(grep -c "^+$" "$FILTERED_FASTQ")
AVG_READ_LEN=$(awk '(NR-2)%4==0 { total += length($0); count++ } END { if (count>0) print total/count; else print 0 }' "$FILTERED_FASTQ")

echo "Read count: $NUM_READS"
echo "Average read length: $AVG_READ_LEN"

# # Step 3: Filter reads if too short
# if [ "$(printf "%.0f" "$AVG_READ_LEN")" -lt 50 ]; then
#     echo "Filtering out reads shorter than 50bp..."
#     seqtk seq -L 50 "$FILTERED_FASTQ" > "filtered_high_quality.fastq"
#     FILTERED_FASTQ="filtered_high_quality.fastq"
#     AVG_READ_LEN=$(awk '(NR-2)%4==0 { total += length($0); count++ } END { if (count>0) print total/count; else print 0 }' "$FILTERED_FASTQ")
#     echo "New average read length: $AVG_READ_LEN"
# fi

# Step 4: Check read quality
if [ "$NUM_READS" -lt 500 ]; then
    echo "Error: Not enough reads for assembly. Skipping." >&2
    exit 1
fi

# Step 5: Run MEGAHIT assembly
echo "Running MEGAHIT..."
rm -rf "$ASSEMBLY_DIR"
megahit -r "$FILTERED_FASTQ" -o "$ASSEMBLY_DIR" --min-contig-len 200

# Step 6: Index reference genome if needed
if [ ! -f "$REFERENCE.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$REFERENCE"
fi

# Step 7: Align contigs back to reference
if [ -f "$CONTIGS_FASTA" ]; then
    echo "Aligning assembled contigs back to reference..."
    bwa mem "$REFERENCE" "$CONTIGS_FASTA" > "$CONTIG_ALIGN_SAM"
else
    echo "Error: MEGAHIT contigs not found at $CONTIGS_FASTA" >&2
    echo "Check MEGAHIT log: $ASSEMBLY_DIR/log" >&2
    exit 1
fi

echo "Done. Contigs assembled and aligned."
echo "Check $CONTIGS_FASTA for assembled sequences"
echo "Check $CONTIG_ALIGN_SAM for mapping results"

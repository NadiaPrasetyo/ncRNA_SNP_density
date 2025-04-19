#!/bin/bash

# Set paths (customize these if needed)
FILTERED_BAM="output/SRR015379_2.filtered.bam"
FILTERED_FASTQ="output/filtered_reads.fastq"
ASSEMBLY_DIR="output/spades_filtered_assembly"
REFERENCE="GRCh38.primary_assembly.genome.fa"
CONTIGS_FASTA="$ASSEMBLY_DIR/contigs.fasta"
CONTIG_ALIGN_SAM="output/contigs_vs_reference.sam"

# Step 1: Convert filtered BAM to FASTQ
# This extracts the actual read sequences from the filtered alignments.
echo "Converting filtered BAM to FASTQ..."
samtools fastq "$FILTERED_BAM" -o "$FILTERED_FASTQ"

# Step 2: Run SPAdes de novo assembly on the extracted reads
# Use lower threshold for coverage and skip error correction to be more permissive.
echo "Running SPAdes for local assembly with lower thresholds..."
spades.py -s "$FILTERED_FASTQ" -o "$ASSEMBLY_DIR" --only-assembler --cov-cutoff 0 --k 21,33,55

# Step 3: Index the reference genome (if not already indexed)
# Required for aligning assembled contigs back to the genome.
if [ ! -f "$REFERENCE.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$REFERENCE"
fi

# Step 4: Align contigs back to the reference genome
# This helps detect if the assembled sequences exist in multiple places.
echo "Aligning assembled contigs back to reference..."
bwa mem "$REFERENCE" "$CONTIGS_FASTA" > "$CONTIG_ALIGN_SAM"

# Optional: Convert the SAM to BAM and sort/index if you want to visualize in IGV
# samtools view -b "$CONTIG_ALIGN_SAM" | samtools sort -o output/contigs_vs_reference.sorted.bam
# samtools index output/contigs_vs_reference.sorted.bam

echo "Done. Contigs assembled and aligned."
echo "Check $CONTIGS_FASTA for assembled sequences"
echo "Check $CONTIG_ALIGN_SAM for where they map in the genome"

# Optionally: Use BLAST for more flexible multi-location matching (commented out by default)
# echo "Running BLAST (requires pre-built BLAST DB from reference genome)..."
# blastn -query "$CONTIGS_FASTA" -db GRCh38.fa -out output/blast_results.txt -outfmt 6

def get_chromosome_locations(genome_fasta):
    """
    Get the start and end positions of each chromosome in the genome FASTA file.
    Returns a dictionary mapping chromosome names to their start and end positions.
    """
    chromosome_locations = {}
    current_chrom = None
    start_pos = 0
    current_pos = 0

    with open(genome_fasta, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                if current_chrom is not None:
                    # Record the end position of the previous chromosome
                    chromosome_locations[current_chrom] = (start_pos, current_pos - 1)
                current_chrom = line[1:]  # Get chromosome name
                start_pos = current_pos + 1  # Start position of the new chromosome
            else:
                current_pos += len(line)  # Increment position by line length

        # Add the last chromosome
        if current_chrom is not None:
            chromosome_locations[current_chrom] = (start_pos, current_pos)

    # Print all chromosome ranges
    print("Chromosome ranges:")
    for chrom, (start, end) in chromosome_locations.items():
        print(f"{chrom}: {start} - {end}")

    return chromosome_locations


def calculate_flanking_n_counts(genome_fasta):
    """
    Calculate the number of leading 'N' nucleotides for each chromosome in the genome.
    Returns a dictionary mapping chromosome names to their N-flank counts.
    """
    flanking_n_counts = {}
    current_chrom = None

    with open(genome_fasta, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_chrom = line[1:]  # Get chromosome name
                flanking_n_counts[current_chrom] = 0
            elif current_chrom:
                # Count leading Ns only for the first sequence line of the chromosome
                for base in line.upper():
                    if base == 'N':
                        flanking_n_counts[current_chrom] += 1
                    else:
                        break

    return flanking_n_counts


def adjust_sam_positions(input_sam, output_sam, chromosome_locations, flanking_n_counts):
    """
    Adjust SAM file positions by mapping reads to chromosomes based on absolute positions,
    considering flanking N counts.
    Writes the adjusted SAM file to the output.
    """
    with open(input_sam, 'r') as infile, open(output_sam, 'w') as outfile:
        for line in infile:
            if line.startswith('@'):  # SAM header lines, write them directly
                outfile.write(line)
                continue

            fields = line.strip().split('\t')
            read_pos = int(fields[4])  # Absolute position from SAM
                        
            mapped_chrom = None
            chrom_pos = None
            
            # Find the chromosome with the smallest difference in position
            for chrom, (start, end) in chromosome_locations.items():
                chrom_name = chrom.split("_")[0].replace("chr", "")
                if chrom_name.isdigit():
                    number = int(chrom_name)
                    if number >1:
                        read_pos = read_pos - number * 2000000
                    
                    if read_pos <= 0:
                        read_pos = read_pos + number * 2000000
                        

                
                if start <= read_pos <= end:
                    mapped_chrom = chrom
                    chrom_pos = read_pos - 1  # Position relative to chromosome
                    # Debug: Successfully mapped the position
                    break
                else:
                    if number>1:
                        read_pos = read_pos + number * 2000000

            if mapped_chrom and mapped_chrom in flanking_n_counts:
                # Adjust position based on flanking N counts
                adjusted_pos = chrom_pos - flanking_n_counts[mapped_chrom]
                if adjusted_pos < 1:  # Ensure positions don't go below 1
                    adjusted_pos = 1
                fields[3] = str(adjusted_pos)  # Update position field
                fields[2] = mapped_chrom  # Update chromosome name
            else:
                print(f"Warning: Position {read_pos} could not be mapped. Line skipped.")
                continue

            outfile.write('\t'.join(fields) + '\n')



def main():
    genome_fasta = "data/animals/hg38.fa"  # Input reference genome file (FASTA)
    input_sam = "data/aligned_genome_reads.sam"  # Input SAM file
    output_sam = "data/adjusted_genome_alignment.sam"  # Output SAM file with adjusted positions

    print("Determining chromosome locations...")
    chromosome_locations = get_chromosome_locations(genome_fasta)

    print("Calculating flanking N counts...")
    flanking_n_counts = calculate_flanking_n_counts(genome_fasta)

    print("Adjusting SAM file positions...")
    adjust_sam_positions(input_sam, output_sam, chromosome_locations, flanking_n_counts)
    print(f"Adjusted SAM file written to {output_sam}")


if __name__ == "__main__":
    main()

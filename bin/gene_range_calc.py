import statistics

def parse_sam_file(sam_file):
    """
    Parses a SAM file to extract positions for each gene.

    Parameters:
        sam_file (str): Path to the input SAM file.

    Returns:
        dict: A dictionary where keys are gene names, and values are lists of positions.
    """
    gene_positions = {}

    with open(sam_file, 'r') as file:
        for line in file:
            if line.startswith("@"):  # Skip SAM headers
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue  # Skip malformed lines
            
            gene = fields[1]
            position = fields[4]

            try:
                position = int(position)  # Convert position to integer
                if gene not in gene_positions:
                    gene_positions[gene] = []
                gene_positions[gene].append(position)
            except ValueError:
                continue  # Skip invalid positions
    
    return gene_positions

def calculate_statistics(gene_positions):
    """
    Calculates the range, mean, and standard deviation of positions for each gene.

    Parameters:
        gene_positions (dict): A dictionary with gene names as keys and positions as values.

    Returns:
        dict: A dictionary containing range, mean, and standard deviation for each gene.
    """
    statistics_summary = {}

    for gene, positions in gene_positions.items():
        if len(positions) < 2:
            statistics_summary[gene] = {"range": 0, "mean": positions[0] if positions else 0, "stdev": 0}
            continue

        max_pos = max(positions)
        min_pos = min(positions)
        pos_range = max_pos - min_pos
        mean_pos = statistics.mean(positions)
        stdev_pos = statistics.stdev(positions)

        statistics_summary[gene] = {
            "range": pos_range,
            "mean": mean_pos,
            "stdev": stdev_pos
        }
    
    return statistics_summary

def main():
    # Hardcoded file name
    sam_file = "data/aligned_genome_reads.sam"  # Replace 'input.sam' with your actual SAM file path
    
    print(f"Processing file: {sam_file}")
    gene_positions = parse_sam_file(sam_file)
    stats = calculate_statistics(gene_positions)

    print(f"{'Gene':<20} {'Range':<10} {'Mean':<15} {'Standard Deviation':<20}")
    print("=" * 65)
    for gene, stat in stats.items():
        print(f"{gene:<20} {stat['range']:<10} {stat['mean']:<15.2f} {stat['stdev']:<20.2f}")

if __name__ == "__main__":
    main()

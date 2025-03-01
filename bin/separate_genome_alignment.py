import os
import re
from collections import defaultdict

def separate_sam_by_genome_part(input_sam, output_dir):
    """
    Separates alignments in a SAM file into separate files based on the genome part of the filename.
    Also counts genes that are not found in all genomes and collects genome appearances.

    Args:
        input_sam (str): Path to the input SAM file.
        output_dir (str): Directory where the separated files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store open file handles for each genome
    file_handles = {}

    # Regular expression to extract genome from the filename
    genome_pattern = re.compile(r"-(\w+)\.sam$")

    # Dictionary to track genes and the genomes they appear in
    gene_presence = defaultdict(set)
    all_genomes = set()

    try:
        with open(input_sam, 'r') as infile:
            for line in infile:
                # Skip headers (lines starting with '@')
                if line.startswith('@'):
                    continue
                
                # Split the line by tab to extract the filename and the rest of the line
                fields = line.strip().split('\t')
                full_filename = fields[0]  # First column is the filename (assumed to contain gene info)

                # Extract genome part using regex
                match = genome_pattern.search(full_filename)
                if not match:
                    continue  # Skip lines with malformed filenames
                
                genome_part = match.group(1)  # Extracted genome name
                all_genomes.add(genome_part)

                # Extract gene name (assuming it's in the first field)
                gene_name = full_filename.split('-')[0]  # Adjust if necessary

                # Track gene presence in genomes
                gene_presence[gene_name].add(genome_part)

                # Remove the filename column to preserve original SAM format
                sam_line = '\t'.join(fields[1:])
                
                # Open a file handle for this genome part if not already opened
                if genome_part not in file_handles:
                    output_path = os.path.join(output_dir, f"{genome_part}.sam")
                    file_handles[genome_part] = open(output_path, 'w')
                
                # Write the SAM line to the corresponding genome file
                file_handles[genome_part].write(sam_line + '\n')

    except Exception as e:
        print(f"Error: {e}")
    
    finally:
        # Close all open file handles
        for handle in file_handles.values():
            handle.close()
        print("Separation complete. All files have been saved.")

    # Identify genes not found in all genomes
    missing_genes = {
        gene: genomes for gene, genomes in gene_presence.items()
        if genomes != all_genomes
    }

    # Print missing gene information
    print("\nGenes not found in all genomes:")
    for gene, genomes in missing_genes.items():
        missing_in = all_genomes - genomes
        print(f"Gene: {gene} | Present in: {sorted(genomes)} | Missing in: {sorted(missing_in)}")

if __name__ == "__main__":
    # Input SAM file
    input_sam = "data/aligned_genome_reads.sam"  # Replace with your input SAM file path
    
    # Output directory to save separated SAM files
    output_dir = "data/separated_genome_sam"  # Replace with your output directory
    
    # Run the function
    separate_sam_by_genome_part(input_sam, output_dir)

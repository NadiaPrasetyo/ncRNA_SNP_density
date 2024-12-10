import pandas as pd
import os
import psutil

def log_memory_usage(stage):
    """Logs the current memory usage for diagnostics."""
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / (1024 * 1024)  # Convert bytes to MB
    print(f"[{stage}] Memory Usage: {memory_mb:.2f} MB")

def estimate_chunk_size(input_vcf, columns, memory_limit_mb=4000):
    """
    Estimate the maximum chunk size to fit within the given memory limit.
    """
    sample_size = 1000
    sample_data = pd.read_csv(input_vcf, sep='\t', comment='#', names=columns, skiprows=100, nrows=sample_size)
    memory_per_row = sample_data.memory_usage(index=True, deep=True).sum() / sample_size
    max_chunk_size = int((memory_limit_mb * 1024 * 1024) / memory_per_row)
    print(f"Estimated maximum chunk size: {max_chunk_size} rows (approx. {memory_limit_mb} MB)")
    return max_chunk_size

def filter_vcf_by_genome(input_vcf, output_dir, memory_limit_mb=4000):
    """
    Filters a VCF file to create individual VCF files for each genome,
    optimized for limited memory.
    
    Args:
        input_vcf (str): Path to the input VCF file.
        output_dir (str): Directory where the filtered VCF files will be saved.
        memory_limit_mb (int): Approximate memory limit for processing in MB.
    """
    log_memory_usage("Start")
    
    # Read header lines
    header_lines = []
    with open(input_vcf, 'r') as file:
        for line in file:
            if line.startswith("##"):
                # skip the line
                continue
            elif line.startswith("#"):
                header_line = line.strip()  # This contains the column headers
                break

    # Parse column names
    column_names = header_line.lstrip("#").split("\t")
    genomes = column_names[10:]  # Extract genomes columns after FORMAT

    # Estimate chunk size based on memory limit
    chunk_size = estimate_chunk_size(input_vcf, column_names, memory_limit_mb)
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    # Process the VCF file in chunks
    data_chunks = pd.read_csv(
        input_vcf, sep='\t', comment='#', names=column_names, skiprows=len(header_lines),
        chunksize=chunk_size
    )

    for chunk_idx, chunk in enumerate(data_chunks):
        log_memory_usage(f"Processing Chunk {chunk_idx}")
        chunk.columns = chunk.columns.str.strip()  # Ensure column names are clean

        for genome in genomes:
            print(f"Processing genome: {genome}, Chunk {chunk_idx}")
            
            # Filter variations
            chunk[genome] = chunk[genome].astype(str)
            genome_variations = chunk[chunk[genome].apply(
                lambda x: any(
                    allele not in {'0', '.', '|'} for allele in x.split('|')
                )
            )]

            if genome_variations.empty:
                continue

            # Select necessary columns
            trimmed_columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', genome]
            genome_variations = genome_variations[trimmed_columns]
            
            # Adjust ALT values on-demand
            def adjust_alt(row):
                alt_alleles = row['ALT'].split(',')
                genotype = row[genome]
                relevant_indices = [
                    int(allele) - 1
                    for allele in genotype.replace('|', '/').split('/')
                    if allele.isdigit() and int(allele) > 0
                ]
                filtered_alt = [alt_alleles[i] for i in relevant_indices]
                return ','.join(filtered_alt) if filtered_alt else '.'

            genome_variations['ALT'] = genome_variations.apply(adjust_alt, axis=1)

            # Append results to output file
            output_file = f"{output_dir}/{genome}.vcf"
            with open(output_file, 'a') as out_vcf:
                if chunk_idx == 0:  # Write headers only for the first chunk
                    genome_header = header_lines[:]
                    genome_header.append(f"##Filtered for genome {genome}\n")
                    out_vcf.writelines(genome_header)
                    out_vcf.write("#" + "\t".join(genome_variations.columns) + "\n")
                genome_variations.to_csv(out_vcf, sep='\t', index=False, header=False)
        
        log_memory_usage(f"After Processing Chunk {chunk_idx}")

    log_memory_usage("End")

# Usage
filter_vcf_by_genome("data/decomposed.vcf", "data/VCF", memory_limit_mb=4000)

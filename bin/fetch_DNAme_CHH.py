import os
import csv
import re
import logging

# Set up logging for detailed diagnostics
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Function to process .bed files in chunks
def process_bed_files(bed_folder, gene_csv, output_csv, chunk_size=1000):
    logging.info("Reading gene ranges from CSV...")

    # Create a dictionary for genes by chromosome
    gene_ranges_by_chromosome = {}
    with open(gene_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for i, row in enumerate(reader):
            if i >= 250:
                break  # Stop reading after 250 lines
            
            chromosome = row['Chromosome']
            if chromosome not in gene_ranges_by_chromosome:
                gene_ranges_by_chromosome[chromosome] = []
            gene_ranges_by_chromosome[chromosome].append({
                'Start': int(row['Start']),
                'End': int(row['End']),
                'GeneName': row['GeneName'],
                'Strand': row.get('Strand', '.'),
                'Length': row['Length'],
                'Median_CG_Content': row['Median_CG_Content']
            })
    
    logging.info(f"Loaded gene ranges for {len(gene_ranges_by_chromosome)} chromosomes.")

    # Prepare output CSV
    output_fields = ['Chromosome', 'Strand', 'Start', 'End', 'Methylation_Percentage', 'GeneName', 'Tissue', 'Length', 'Median_CG_Content']

    with open(output_csv, 'a', newline='') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=output_fields)
        writer.writeheader()  # Writing the header only once

        # Process each .bed file
        bed_files = sorted([f for f in os.listdir(bed_folder) if f.endswith('.bed')])
        for bed_file in bed_files:
            match = re.search(r'_([a-zA-Z0-9]+)_CHH\.bed$', bed_file)
            tissue_name = match.group(1) if match else "unknown"

            # skip if the files are epithelium, brain, or embryo
            if tissue_name in ['epithelium', 'brain', 'embryo']:
                continue


            bed_file_path = os.path.join(bed_folder, bed_file)
            logging.info(f"Processing {bed_file} (tissue: {tissue_name})...")

            processed_lines = 0

            try:
                with open(bed_file_path, 'r') as bed:
                    chunk = []
                    for line in bed:
                        processed_lines += 1
                        if len(chunk) >= chunk_size:
                            process_chunk(chunk, gene_ranges_by_chromosome, writer, tissue_name)
                            chunk.clear()
                            logging.info(f"Processed {processed_lines} lines so far in {bed_file}...")

                        fields = line.strip().split('\t')
                        bed_chrom = fields[0]
                        bed_start = int(fields[1])
                        bed_end = int(fields[2])
                        strand = fields[5]
                        methylation_percentage = float(fields[10])

                        chunk.append((bed_chrom, bed_start, bed_end, strand, methylation_percentage))

                    if chunk:
                        process_chunk(chunk, gene_ranges_by_chromosome, writer, tissue_name)
                        logging.info(f"Finished processing {bed_file}.")
            except Exception as e:
                logging.error(f"Error processing {bed_file}: {e}")

    logging.info(f"Methylation data has been written to {output_csv}.")

# Function to process a chunk of lines
def process_chunk(chunk, gene_ranges_by_chromosome, writer, tissue_name):
    overlaps_found = 0
    skipped_lines = 0

    for line in chunk:
        bed_chrom, bed_start, bed_end, strand, methylation_percentage = line

        if bed_chrom not in gene_ranges_by_chromosome:
            skipped_lines += 1
            continue

        for gene in gene_ranges_by_chromosome[bed_chrom]:
            if bed_start <= gene['End'] and bed_end >= gene['Start']:
                overlaps_found += 1
                writer.writerow({
                    'Chromosome': bed_chrom,
                    'Strand': strand,
                    'Start': max(bed_start, gene['Start']),
                    'End': min(bed_end, gene['End']),
                    'Methylation_Percentage': methylation_percentage,
                    'GeneName': gene['GeneName'],
                    'Tissue': tissue_name,
                    'Length': gene['Length'],
                    'Median_CG_Content': gene['Median_CG_Content']
                })
    logging.info(f"Processed chunk with {overlaps_found} overlaps and {skipped_lines} skipped lines.")

# Paths to input and output files
bed_folder_path = "data/"
gene_csv_path = "data/SNP_RNA_GC.csv"
output_csv_path = "data/CHH_methylation_data_temp_temp.csv"

# Run the processing function
process_bed_files(bed_folder_path, gene_csv_path, output_csv_path)

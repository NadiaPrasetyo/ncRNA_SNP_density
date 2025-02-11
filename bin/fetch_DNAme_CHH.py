import os
import csv
import re

# Function to read .bed files and process gene ranges
def process_bed_files(bed_folder, gene_csv, output_csv):
    print("Reading gene ranges from CSV...")

    # Read the gene ranges from the input CSV and organize them by chromosome
    gene_ranges_by_chromosome = {}
    with open(gene_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            chromosome = row['Chromosome']
            if chromosome not in gene_ranges_by_chromosome:
                gene_ranges_by_chromosome[chromosome] = []
            gene_ranges_by_chromosome[chromosome].append({
                'Start': int(row['Start']),
                'End': int(row['End']),
                'GeneName': row['GeneName'],
                'Strand': row.get('Strand', '.'),  # Default strand to '.' if not provided
                'Length': row['Length'],
                'Median_CG_Content': row['Median_CG_Content']
            })

    chromosomes_to_process = set(gene_ranges_by_chromosome.keys())
    print(f"Loaded gene ranges for {len(chromosomes_to_process)} chromosomes.")

    # Prepare the output CSV
    output_fields = ['Chromosome', 'Strand', 'Start', 'End', 'Methylation_Percentage', 'GeneName', 'Tissue', 'Length', 'Median_CG_Content']
    with open(output_csv, 'w', newline='') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=output_fields)
        writer.writeheader()

        # Process each .bed file
        bed_files = [f for f in os.listdir(bed_folder) if f.endswith('.bed')]
        print(f"Found {len(bed_files)} .bed files to process.")

        for bed_file in bed_files:
            # Extract tissue name from the file name using regex
            match = re.search(r'_([a-zA-Z0-9]+)_CHH\.bed$', bed_file)
            tissue_name = match.group(1) if match else "unknown"

            bed_file_path = os.path.join(bed_folder, bed_file)

            print(f"Processing {bed_file} (tissue: {tissue_name})...")
            processed_lines = 0
            overlaps_found = 0
            skipped_lines = 0

            try:
                with open(bed_file_path, 'r') as bed:
                    for line in bed:
                        processed_lines += 1

                        # Process every 1000 lines to monitor progress
                        if processed_lines % 1000 == 0:
                            print(f"Processed {processed_lines} lines so far in {bed_file}...")

                        fields = line.strip().split('\t')
                        bed_chrom = fields[0]
                        bed_start = int(fields[1])
                        bed_end = int(fields[2])
                        strand = fields[5]
                        methylation_percentage = float(fields[10])

                        # Skip lines for chromosomes not in the CSV
                        if bed_chrom not in chromosomes_to_process:
                            skipped_lines += 1
                            continue

                        # Check if this region overlaps with any gene range for this chromosome
                        for gene in gene_ranges_by_chromosome[bed_chrom]:
                            if bed_start <= gene['End'] and bed_end >= gene['Start']:
                                overlaps_found += 1
                                # Collect overlapping methylation data
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

                print(f"Processed {processed_lines} lines from {bed_file}. Found {overlaps_found} overlaps. Skipped {skipped_lines} lines.")
            except Exception as e:
                print(f"Error processing {bed_file}: {e}")

    print(f"Methylation data has been written to {output_csv}.")

# Paths to input and output files
bed_folder_path = "data/"
gene_csv_path = "data/SNP_RNA_GC.csv"
output_csv_path = "data/CHH_methylation_data.csv"

# Call the function to process the files
process_bed_files(bed_folder_path, gene_csv_path, output_csv_path)


#order: chromosome, start, end, name of item, score, strand, colour, coverage/number of reads, percentage of reads that show methylation at this position in the genome
#Processed 55883000 lines so far in ENCFF839MWC_brain_CHH.bed...
# Processed 55884000 lines so far in ENCFF839MWC_brain_CHH.bed...
# Processed 55885000 lines so far in ENCFF839MWC_brain_CHH.bed...
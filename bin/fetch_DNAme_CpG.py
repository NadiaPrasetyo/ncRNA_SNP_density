import os
import csv
import pyBigWig
import traceback

# Directory where your BigBed files are stored
data_directory = 'data/'

# Path to your CSV file containing gene regions
input_csv_file = 'data/SNP_RNA_GC.csv'

# Output CSV file path
output_csv_file = 'output_methylation_data.csv'

# Function to process the BigBed files and fetch methylation data for the given gene regions
def process_bigbed_file_for_genes(bigbed_filename, genes, output_writer):
    try:
        print(f"\nProcessing file: {bigbed_filename}")
        
        # Extract tissue from the filename (assumes format: ENCFFXXXXXX_tissue_CpG.bigBed)
        tissue = os.path.basename(bigbed_filename).split('_')[1]
        
        # Open the BigBed file with pyBigWig
        bw = pyBigWig.open(bigbed_filename)
        
        if bw is None:
            print(f"Failed to open {bigbed_filename}: File might not be a valid BigBed file.")
            return
        
        # Iterate through each gene
        for gene in genes:
            chrom = gene['Chromosome']
            gene_start = gene['Start']
            gene_end = gene['End']
            gene_name = gene['GeneName']

            print(f"  Querying region: {chrom}:{gene_start}-{gene_end} for gene: {gene_name}")
            
            # Fetch the entries from the BigBed file that overlap the gene's region
            entries = bw.entries(chrom, gene_start, gene_end)
            
            if entries is None:
                print(f"  No entries found for region: {chrom}:{gene_start}-{gene_end}")
                continue
            
            for start, end, score in entries:
                # Process the tab-separated data in score
                tab_data = score.split('\t')
                
                # Check if the tab data has the expected number of fields
                if len(tab_data) < 8:
                    print(f"  Malformed data in entry: {tab_data}")
                    continue
                
                # Extract the relevant fields from the tab-separated data
                strand = tab_data[2]
                methylation_percentage = tab_data[7]  # Percentage of read showing methylation

                # Write the relevant information to the CSV
                output_writer.writerow([chrom, strand, start, end, methylation_percentage, gene_name, tissue])

        # Close the BigBed file
        bw.close()

    except Exception as e:
        print(f"Error processing file {bigbed_filename}: {e}")
        print("Traceback:")
        traceback.print_exc()

# Function to read the CSV file containing gene regions
def read_gene_regions(csv_filename):
    genes = []
    try:
        with open(csv_filename, mode='r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                # Extract the necessary fields
                gene = {
                    'Chromosome': row['Chromosome'],
                    'Start': int(row['Start']),
                    'End': int(row['End']),
                    'GeneName': row['GeneName']  # Add GeneName for output
                }
                genes.append(gene)
    except Exception as e:
        print(f"Error reading CSV file {csv_filename}: {e}")
        print("Traceback:")
        traceback.print_exc()
    return genes

# Main execution
def main():
    # Read gene regions from the input CSV file
    genes = read_gene_regions(input_csv_file)
    
    if not genes:
        print("No genes loaded from the CSV file. Exiting.")
        return
    
    # Open the output CSV file for writing
    with open(output_csv_file, mode='w', newline='') as output_file:
        output_writer = csv.writer(output_file)
        # Write the header row
        output_writer.writerow(['Chromosome', 'Strand', 'Start', 'End', 'Methylation_Percentage', 'GeneName', 'Tissue'])
        
        # Loop through all the BigBed files in the data directory
        for filename in os.listdir(data_directory):
            if filename.endswith('.bigBed'):  # Check if it's a BigBed file
                process_bigbed_file_for_genes(os.path.join(data_directory, filename), genes, output_writer)

# Run the script
if __name__ == "__main__":
    main()

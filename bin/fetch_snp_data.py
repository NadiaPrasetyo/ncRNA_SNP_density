import pyBigWig
import pandas as pd

# Path to your BigBed file
bigbed_file = "data/dbSnp155.bb"

# Path to your CSV file containing genes with chromosomal locations
csv_file = "data/SNP-densities-and-RNA.csv"  # CSV should have 'gene_name', 'chrom', 'start', 'end' columns

# Output text file where data will be saved
output_file = "data/SNP_data.txt"


# Open the output file in write mode to store the results
with open(output_file, 'w') as f:
    # Path to BigBed file
    bb = pyBigWig.open(bigbed_file)

    try:
        if bb.isBigBed():
            print("Opened BigBed file successfully.")
            print(f"Header Info: {bb.header()}\n")
            f.write("SNP details within gene regions:\n\n")  # Start output file with a header for clarity

            # Load gene data from CSV file
            gene_data = pd.read_csv(csv_file)

            # Debugging: Print column names in the CSV file (for diagnostics)
            print(f"CSV Columns: {gene_data.columns.tolist()}")

            # Iterate through the CSV file (one gene at a time)
            for _, gene_row in gene_data.iterrows():
                chrom = gene_row['Chromosome']   # Chromosome name (e.g., 'chr1')
                start = gene_row['Start']        # Gene start position
                end = gene_row['End']           # Gene end position
                gene_name = gene_row['GeneName'] # Gene name (e.g., 'RNU5E-1')

                # Print gene, chromosome, start, end for diagnostics
                print(f"Processing gene: {gene_name} on {chrom} from {start} to {end}")

                # Fetch intervals from the BigBed file for this specific region
                intervals = bb.entries(chrom, start, end)

                # List to store SNP names for the current gene
                snp_names = []

                if intervals:
                    for interval in intervals:
                        # The interval is a tuple with the format: (chromosome, start_position, data)
                        snp_details = interval[2].split("\t")  # Split the details string by tab

                        try:
                            chromStart = interval[0]  # Start position
                            chromEnd = interval[1]  # End position (next base after start for single-base variants)
                            name = snp_details[0]  # SNP name (dbSNP ID)
                            ref = snp_details[1]  # Reference allele
                            alts = snp_details[3].split(",")[:int(snp_details[2])]  # Alternate alleles
                            variation_class = snp_details[10]  # Variation class (e.g., SNV)
                            ucsc_notes = snp_details[11]  # UCSC notes

                            # Collect SNP names for diagnostics
                            snp_names.append(name)

                            # Write the required information to the output file
                            f.write(f"Gene: {gene_name}\n")
                            f.write(f"Chromosome: {chrom}\n")
                            f.write(f"Start: {chromStart}, End: {chromEnd}\n")
                            f.write(f"SNP Name: {name}\n")
                            f.write(f"Reference Allele: {ref}\n")
                            f.write(f"Alternate Alleles: {', '.join(alts)}\n")
                            f.write(f"Variation Class: {variation_class}\n")
                            f.write(f"UCSC Notes: {ucsc_notes}\n")
                            f.write("-" * 40 + "\n")

                        except IndexError:
                            print(f"Error processing SNP details for: {interval}")

                # After processing the SNPs for the current gene, print the diagnostics
                print(f"Gene: {gene_name} - Chromosome: {chrom} | Start: {start}, End: {end}")
                print(f"SNPs found: {', '.join(snp_names)}\n")

        else:
            print("The file is not a BigBed file.")

    except Exception as e:
        print(f"An error occurred: {e}")

    finally:
        # Close the BigBed file
        bb.close()

    print("Processing completed.")
    print(f"SNP data has been written to {output_file}.")
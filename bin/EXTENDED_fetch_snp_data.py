import pyBigWig
import pandas as pd

# Path to your BigBed file
bigbed_file = "data/dbSnp155.bb"

# Path to your CSV file containing genes with chromosomal locations
csv_file = "data/SNP-densities-and-RNA.csv"  # CSV should have 'gene_name', 'chrom', 'start', 'end' columns

# Output text file where data will be saved
output_file = "data/EXTENDED_SNP_data.txt"

# Open the output file in write mode to store the results
with open(output_file, 'w') as f:
    # Path to BigBed file
    bb = pyBigWig.open(bigbed_file)

    try:
        if bb.isBigBed():
            print("Opened BigBed file successfully.")
            print(f"Header Info: {bb.header()}\n")
            
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

                # Calculate the expanded boundaries
                gene_length = end - start
                flank_length = int(4.5 * gene_length)
                expanded_start = max(0, start - flank_length)  # Avoid negative start positions
                expanded_end = end + flank_length

                # Print gene, chromosome, and expanded boundaries for diagnostics
                print(f"Processing gene: {gene_name} on {chrom}")
                print(f"Original: {start}-{end} | Expanded: {expanded_start}-{expanded_end}")

                # Fetch intervals from the BigBed file for the expanded region
                intervals = bb.entries(chrom, expanded_start, expanded_end)

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

                            # Check if the SNP is in the gene or flank region
                            in_gene = chromStart >= start and chromEnd <= end
                            in_flank = chromStart >= expanded_start and chromEnd <= expanded_end

                            if in_gene and in_flank:
                                # SNP is in both the gene and flank region, label only with gene_name
                                region_type = gene_name
                            elif in_gene:
                                # SNP is inside the gene
                                region_type = gene_name
                            elif in_flank:
                                # SNP is in the flank region
                                region_type = f"{gene_name}_flank"

                            # Collect SNP names for diagnostics
                            snp_names.append(name)

                            # Write the required information to the output file
                            f.write(f"Gene: {gene_name}\n")
                            f.write(f"Chromosome: {chrom}\n")
                            f.write(f"Expanded Start: {expanded_start}, Expanded End: {expanded_end}\n")
                            f.write(f"Original Start: {start}, Original End: {end}\n")
                            f.write(f"SNP Name: {name}\n")
                            f.write(f"Reference Allele: {ref}\n")
                            f.write(f"Alternate Alleles: {', '.join(alts)}\n")
                            f.write(f"Variation Class: {variation_class}\n")
                            f.write(f"UCSC Notes: {ucsc_notes}\n")
                            f.write(f"Region Type: {region_type}\n")  # Add region type (gene or gene_flank)
                            f.write("-" * 40 + "\n")

                        except IndexError:
                            print(f"Error processing SNP details for: {interval}")

                # After processing the SNPs for the current gene, print the diagnostics
                print(f"Gene: {gene_name} - Chromosome: {chrom} | Expanded Start: {expanded_start}, Expanded End: {expanded_end}")
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

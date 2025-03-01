import pyBigWig
import pandas as pd

# File paths
bigbed_file = "data/dbSnp155.bb"  # Path to your BigBed file
csv_file = "data/SNP-densities-and-RNA.csv"  # Path to CSV containing gene data
output_file = "data/SNP_data.txt"  # Output file for SNP data

def process_gene_snp_data(bigbed_file, csv_file, output_file):
    """
    Processes SNP data for genes listed in a CSV file using a BigBed file.

    Parameters:
    bigbed_file (str): Path to the BigBed file.
    csv_file (str): Path to the CSV file containing gene information.
    output_file (str): Path to the output file where results will be saved.

    Outputs:
    Writes SNP data to the specified output file.
    """
    try:
        # Open the BigBed file
        bb = pyBigWig.open(bigbed_file)
        if not bb.isBigBed():
            raise ValueError("The specified file is not a valid BigBed file.")

        print("Opened BigBed file successfully.")
        print(f"Header Info: {bb.header()}")

        # Load gene data from the CSV file
        gene_data = pd.read_csv(csv_file)

        # Check for required columns
        required_columns = {'Chromosome', 'Start', 'End', 'GeneName'}
        if not required_columns.issubset(gene_data.columns):
            raise KeyError(f"The input CSV file must contain the columns: {', '.join(required_columns)}.")

        print(f"CSV Columns: {gene_data.columns.tolist()}")

        # Open the output file for writing
        with open(output_file, 'w') as f:
            # Process each gene
            for _, gene_row in gene_data.iterrows():
                try:
                    chrom = str(gene_row['Chromosome']).strip()
                    start = int(gene_row['Start'])
                    end = int(gene_row['End'])
                    gene_name = gene_row['GeneName']

                    print(f"Processing gene: {gene_name} on {chrom} from {start} to {end}")

                    # Fetch intervals from the BigBed file
                    intervals = bb.entries(chrom, start, end)
                    snp_names = []

                    if intervals:
                        for interval in intervals:
                            try:
                                # Extract interval details
                                snp_details = interval[2].split("\t")
                                chromStart = interval[0]
                                chromEnd = interval[1]
                                name = snp_details[0]
                                ref = snp_details[1]
                                alts = snp_details[3].split(",")[:int(snp_details[2])]
                                variation_class = snp_details[10]
                                ucsc_notes = snp_details[11]

                                snp_names.append(name)

                                # Write details to the output file
                                f.write(f"Gene: {gene_name}\n")
                                f.write(f"Chromosome: {chrom}\n")
                                f.write(f"Start: {chromStart}, End: {chromEnd}\n")
                                f.write(f"SNP Name: {name}\n")
                                f.write(f"Reference Allele: {ref}\n")
                                f.write(f"Alternate Alleles: {', '.join(alts)}\n")
                                f.write(f"Variation Class: {variation_class}\n")
                                f.write(f"UCSC Notes: {ucsc_notes}\n")
                                f.write("-" * 40 + "\n")
                            except IndexError as ie:
                                print(f"Error processing SNP details for interval: {interval} | Error: {ie}")
                    else:
                        print(f"No SNPs found for gene {gene_name} on {chrom}:{start}-{end}")

                    print(f"SNPs found for {gene_name}: {', '.join(snp_names) if snp_names else 'None'}\n")

                except Exception as gene_error:
                    print(f"Error processing gene {gene_row['GeneName']}: {gene_error}")

        print(f"SNP data has been written to {output_file}.")

    except FileNotFoundError as fnfe:
        print(f"Error: {fnfe}")
    except ValueError as ve:
        print(f"Error: {ve}")
    except KeyError as ke:
        print(f"Error: {ke}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        try:
            # Close the BigBed file if it was opened
            bb.close()
            print("BigBed file closed successfully.")
        except:
            print("Failed to close the BigBed file.")

if __name__ == "__main__":
    process_gene_snp_data(bigbed_file, csv_file, output_file)

import pandas as pd
import pyBigWig
import random

# Input and output file paths
input_csv = "data/SNP-densities-and-RNA.csv"
output_csv = "data/SNP_RNA_GC.csv"
gc_content_bw = "data/hg38.gc5Base.bw"

# Function to calculate median CG content for a genomic region
def calculate_median_gc(chromosome, start, end, bw):
    try:
        values = bw.values(chromosome, start, end, numpy=True)
        # Remove None values (regions with no data)
        values = values[~pd.isnull(values)]
        if len(values) > 0:
            return float(pd.Series(values).median())
        else:
            return None
    except Exception as e:
        print(f"Error processing {chromosome}:{start}-{end}: {e}")
        return None

# Function to calculate the length of a genomic region
def calculate_length(start, end):
    return end - start

def generate_random_regions(length_range, num_regions, bw):
    random_regions = []
    canonical_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']  # Canonical chromosomes (chr1-22, chrX, chrY)    
    bw_chroms = bw.chroms()  # Get available chromosomes in BigWig file
    
    # Filter to include only chromosomes that are canonical and exist in the BigWig file
    valid_chromosomes = [chrom for chrom in canonical_chromosomes if chrom in bw_chroms]

    if not valid_chromosomes:
        print("No valid canonical chromosomes found in the BigWig file.")
        return []

    for i in range(num_regions):
        chrom = random.choice(valid_chromosomes)  # Choose only from valid canonical chromosomes
        chrom_length = bw_chroms[chrom]
        if chrom_length <= 0:  # Skip invalid chromosome lengths
            continue
        length = random.randint(length_range[0], length_range[1])
        max_end = chrom_length - length
        if max_end > 0:
            start = random.randint(0, max_end)
            end = start + length
            median_cg = calculate_median_gc(chrom, start, end, bw)
            random_regions.append({
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Length": length,
                "Median_CG_Content": median_cg,
                "GeneName": f"Random{i+1}",
                "SNP_density": None,
                "rna_type": None
            })
    
    # Ensure we generate exactly the requested number of random regions
    while len(random_regions) < num_regions:
        chrom = random.choice(valid_chromosomes)
        chrom_length = bw_chroms[chrom]
        if chrom_length <= 0:
            continue
        length = random.randint(length_range[0], length_range[1])
        max_end = chrom_length - length
        if max_end > 0:
            start = random.randint(0, max_end)
            end = start + length
            median_cg = calculate_median_gc(chrom, start, end, bw)
            random_regions.append({
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Length": length,
                "Median_CG_Content": median_cg,
                "GeneName": f"Random{len(random_regions)+1}",
                "SNP_density": None,
                "rna_type": None
            })
    
    return random_regions


def main():
    # Read the input CSV
    df = pd.read_csv(input_csv)

    # Open the BigWig file
    bw = pyBigWig.open(gc_content_bw)

    # Add new columns for length and median CG content
    lengths = []
    median_cgs = []

    for index, row in df.iterrows():
        chromosome = row['Chromosome']
        start = int(row['Start'])
        end = int(row['End'])

        # Calculate length
        length = calculate_length(start, end)
        lengths.append(length)

        # Calculate median CG content
        median_cg = calculate_median_gc(chromosome, start, end, bw)
        median_cgs.append(median_cg)

    # Add the new data to the DataFrame
    df['Length'] = lengths
    df['Median_CG_Content'] = median_cgs

    # Generate random regions
    length_range = (min(lengths), max(lengths))
    random_regions = generate_random_regions(length_range, 1000, bw)

    # Append random regions to the DataFrame
    if random_regions:  # Ensure random regions were generated
        random_df = pd.DataFrame(random_regions)
        df = pd.concat([df, random_df], ignore_index=True)

    # Close the BigWig file
    bw.close()

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_csv, index=False)
    print(f"Updated CSV file saved to {output_csv}")

if __name__ == "__main__":
    main()

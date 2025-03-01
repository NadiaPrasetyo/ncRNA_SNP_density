import pyBigWig
import pandas as pd
import os

# Read the SNP and RNA data from the CSV file
csv_file_path = "data/SNP-densities-and-RNA.csv"
df = pd.read_csv(csv_file_path)

# List of BigBed files in the 'data/' folder (only files ending with .bigBed)
bigbed_files = [f for f in os.listdir('data/') if f.endswith(".bigBed")]

# Function to extract data from BigBed file for each gene region
def get_region_data_from_bigbed(bigbed_file, chrom, start, end):
    try:
        # Open the BigBed file
        bb = pyBigWig.open(bigbed_file)
        
        # Fetch the intervals for the given region
        intervals = bb.intervals(chrom, start, end)
        
        # Close the BigBed file after reading
        bb.close()
        
        # Return the intervals (list of tuples: start, end, value)
        return intervals if intervals else []
    except Exception as e:
        print(f"Error fetching data from {bigbed_file} for {chrom}:{start}-{end}: {e}")
        return []

# Iterate through the dataframe and fetch data from the BigBed files
output = []

for _, row in df.iterrows():
    chrom = row['Chromosome']
    start = row['Start']
    end = row['End']
    gene_name = row['GeneName']
    
    print(f"Processing {gene_name} ({chrom}:{start}-{end})...")

    for bigbed_file in bigbed_files:
        print(f"  Fetching region data from {bigbed_file}...")
        
        # Fetch data for the gene region from each BigBed file
        region_data = get_region_data_from_bigbed(os.path.join('data', bigbed_file), chrom, start, end)
        
        if region_data:
            print(f"    Region data for {gene_name} from {bigbed_file}: {region_data[:5]}...")  # Print first few intervals as sample
            
        # Store the result
        output.append({
            'Chromosome': chrom,
            'Start': start,
            'End': end,
            'GeneName': gene_name,
            'BigBedFile': bigbed_file,
            'RegionData': region_data
        })

# Convert the output to a DataFrame
output_df = pd.DataFrame(output)

# Print the output to verify
print("\nProcessed data:\n")
print(output_df.head())  # Print first few rows of the result

# Save the output as a CSV file
output_df.to_csv('data/region_data.csv', index=False)

print("\nRegion data extraction completed and saved to 'data/region_data.csv'.")

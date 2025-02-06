import pyBigWig
import pandas as pd
import os

# Read the SNP and RNA data from the CSV file
csv_file_path = "data/SNP-densities-and-RNA.csv"
df = pd.read_csv(csv_file_path)

# List of bigWig files in the 'data/' folder (only files ending with .bigWig)
bigwig_files = [f for f in os.listdir('data/') if f.endswith(".bigWig")]

# Function to extract methylation data from bigWig file for each gene region
def get_methylation_signal(bigwig_file, chrom, start, end):
    try:
        # Open the bigWig file
        bw = pyBigWig.open(bigwig_file)
        
        # Fetch the values for the given region
        signal = bw.values(chrom, start, end)
        
        # Close the bigWig file after reading
        bw.close()
        
        # Return the signal as a list of values (could be summed or averaged depending on the requirement)
        return signal
    except Exception as e:
        print(f"Error fetching data from {bigwig_file} for {chrom}:{start}-{end}: {e}")
        return []

# Iterate through the dataframe and fetch methylation signals from the bigWig files
output = []

for _, row in df.iterrows():
    chrom = row['Chromosome']
    start = row['Start']
    end = row['End']
    gene_name = row['GeneName']
    
    print(f"Processing {gene_name} ({chrom}:{start}-{end})...")

    for bigwig_file in bigwig_files:
        print(f"  Fetching methylation signal from {bigwig_file}...")
        
        # Fetch methylation signal for the gene region from each bigWig file
        methylation_signal = get_methylation_signal(os.path.join('data', bigwig_file), chrom, start, end)
        
        if methylation_signal:
            print(f"    Methylation signal for {gene_name} from {bigwig_file}: {methylation_signal[:5]}...")  # Print first few values as sample
            
        # Store the result
        output.append({
            'Chromosome': chrom,
            'Start': start,
            'End': end,
            'GeneName': gene_name,
            'BigWigFile': bigwig_file,
            'MethylationSignal': methylation_signal
        })

# Convert the output to a DataFrame
output_df = pd.DataFrame(output)

# Print the output to verify
print("\nProcessed data:\n")
print(output_df.head())  # Print first few rows of the result

# Save the output as a CSV file
output_df.to_csv('data/methylation_signals.csv', index=False)

print("\nMethylation signal extraction completed and saved to 'data/methylation_signals.csv'.")

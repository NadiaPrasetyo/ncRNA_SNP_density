import pandas as pd

# Read the CSV file
input_file = 'data/SNP-densities-and-RNA.csv'  # Replace with your file name
output_file = 'results/UCSC_links.csv'

# Load the data into a pandas DataFrame
df = pd.read_csv(input_file)

# Create a new column for the UCSC link
def generate_ucsc_link(row):
    chromosome = row['Chromosome']
    start = row['Start']
    end = row['End']
    # Construct UCSC link (you can change hg38 to other versions if needed)
    ucsc_link = f'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chromosome}:{start}-{end}'
    return ucsc_link

df['UCSC_Link'] = df.apply(generate_ucsc_link, axis=1)

# Save the updated dataframe to a new CSV
df.to_csv(output_file, index=False)

print(f"UCSC links have been added and saved to {output_file}")

import pandas as pd

# File paths
input_file = 'data/SNP-densities-and-RNA.csv'  # Replace with your file name
output_file = 'results/UCSC_links.csv'

def generate_ucsc_link(row):
    """
    Generates a UCSC Genome Browser link for a given row in the DataFrame.

    Parameters:
    row (pd.Series): A row from the DataFrame containing 'Chromosome', 'Start', and 'End' columns.

    Returns:
    str: A formatted UCSC Genome Browser link.
    """
    try:
        chromosome = str(row['Chromosome']).strip()  # Ensure chromosome is a string
        start = int(row['Start'])  # Ensure start is an integer
        end = int(row['End'])  # Ensure end is an integer
        
        if start > end:
            raise ValueError(f"Start ({start}) is greater than End ({end}) for Chromosome {chromosome}.")
        
        # Construct UCSC link (defaulting to hg38; adjust as needed)
        ucsc_link = f'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={chromosome}:{start}-{end}'
        return ucsc_link
    
    except (ValueError, TypeError) as e:
        print(f"Error generating UCSC link for row {row.name}: {e}")
        return "Invalid Link"

def main():
    """
    Main function to add UCSC Genome Browser links to a CSV and save the updated data.
    """
    try:
        # Load the data into a pandas DataFrame
        df = pd.read_csv(input_file)
        
        # Validate required columns
        required_columns = {'Chromosome', 'Start', 'End'}
        if not required_columns.issubset(df.columns):
            raise KeyError(f"Input file must contain the columns: {', '.join(required_columns)}.")
        
        # Generate UCSC links
        df['UCSC_Link'] = df.apply(generate_ucsc_link, axis=1)

        # Save the updated DataFrame to a new CSV
        df.to_csv(output_file, index=False)
        print(f"UCSC links have been added and saved to {output_file}")
    
    except FileNotFoundError:
        print(f"Error: The input file '{input_file}' was not found.")
    except KeyError as e:
        print(f"Error: {e}")
    except pd.errors.EmptyDataError:
        print(f"Error: The input file '{input_file}' is empty.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()

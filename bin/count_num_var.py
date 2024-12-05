import pandas as pd
from pathlib import Path

def process_gene_variations(csv_file):
    """
    Process a CSV file containing gene variation data and generate a summary of variation counts and types,
    along with gene and total lengths. If any of the required length fields are missing, set the length to "NA".
    
    Args:
        csv_file (str): Path to the input CSV file.
        
    Returns:
        pd.DataFrame: A DataFrame summarizing the variations for each gene.
    """
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file)
    
    # Initialize the results list
    results = []
    
    # Group data by the 'GENE' column
    for gene, gene_data in data.groupby('GENE'):
        # Initialize a dictionary for counting variation types
        variation_counts = {}
        num_variations = 0
        
        # Get the first row of the gene data to extract lengths
        gene_row = gene_data.iloc[0]
        
        # Check if the required fields are missing or empty, and set lengths accordingly
        try:
            original_start = gene_row['ORIGINAL_START']
            original_end = gene_row['ORIGINAL_END']
            extended_start = gene_row['EXTENDED_START']
            extended_end = gene_row['EXTENDED_END']
            
            # Check for missing values and set lengths to 'NA' if any required field is missing
            if pd.isna(original_start) or pd.isna(original_end):
                original_length = 'NA'
            else:
                original_length = original_end - original_start + 1

            if pd.isna(extended_start) or pd.isna(extended_end):
                total_length = 'NA'
            else:
                total_length = extended_end - extended_start + 1

        except KeyError:
            original_length = 'NA'
            total_length = 'NA'

        # Iterate through each row for this gene
        for _, row in gene_data.iterrows():
            # If 'TYPE' is not null or empty, process the variations
            if pd.notnull(row['TYPE']) and row['TYPE'] != "":
                types = row['TYPE'].split(",")  # Split the 'TYPE' by commas
                num_variations += len(types)  # Add the total number of types for this row
                
                # Count occurrences of each variation type
                for type_ in types:
                    type_ = type_.strip()  # Remove any leading/trailing whitespace
                    if type_ in variation_counts:
                        variation_counts[type_] += 1
                    else:
                        variation_counts[type_] = 1
            else:
                # If 'TYPE' is null or empty, treat it as one variation
                num_variations += 1
        
        # Format the variation counts as a string
        formatted_variation_counts = ", ".join([f"{key}({value})" for key, value in variation_counts.items()])
        
        # Compile results for this gene
        results.append({
            'GENE': gene,
            'NUM_VARIATIONS': num_variations,  # The total count of variations (including 1 for null TYPEs)
            'VARIATION_TYPES': formatted_variation_counts,  # List of variation types with counts
            'GENE_LENGTH': original_length,  # Original gene length (or 'NA' if missing)
            'TOTAL_LENGTH': total_length  # Extended total length (or 'NA' if missing)
        })
    
    # Convert results into a DataFrame
    results_df = pd.DataFrame(results)
    
    # Define output path and ensure the directory exists
    output_path = Path('results/EXTENDED_pangenome_variation_summary.csv')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the results to a new CSV file
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to '{output_path}'")
    
    return results_df

# Example usage
if __name__ == "__main__":
    csv_file = 'data/EXTENDED_pangenome_summary.csv'  # Replace with the actual file path
    summary = process_gene_variations(csv_file)
    print(summary)

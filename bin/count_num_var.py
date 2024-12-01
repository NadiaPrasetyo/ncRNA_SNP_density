import pandas as pd
from pathlib import Path

def process_gene_variations(csv_file):
    """
    Process a CSV file containing gene variation data and generate a summary of variation counts and types.
    
    Args:
        csv_file (str): Path to the input CSV file.
        
    Returns:
        pd.DataFrame: A DataFrame summarizing the variations for each gene.
    """
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the gene and type columns
    gene_column = data['GENE']
    
    #count the number of variations for each gene
    results = []
    
    # Read the CSV file into a pandas DataFrame
    data = pd.read_csv(csv_file)
    
    # Extract the gene and type columns
    gene_column = data['GENE']
    
    #count the number of variations for each gene
    results = []
    
    # Group data by the 'GENE' column
    for gene, gene_data in data.groupby('GENE'):
        # Initialize variables for counting
        variation_types = set()
        num_variations = 0
        
        # Iterate through each row for this gene
        for _, row in gene_data.iterrows():
            # If 'TYPE' is not null or empty, process the variations
            if pd.notnull(row['TYPE']) and row['TYPE'] != "":
                types = row['TYPE'].split(",")  # Split the 'TYPE' by commas
                variation_types.update(types)  # Add the types to the set (unique variations)
                num_variations += len(types)  # Add the total number of types for this row
            else:
                # If 'TYPE' is null or empty, add 1 to the count for this position
                num_variations += 1
        
        # Compile results for this gene
        results.append({
            'GENE': gene,
            'NUM_VARIATIONS': num_variations,  # The total count of variations (including 1 for null TYPEs)
            'VARIATION_TYPES': ",".join(variation_types)  # List of unique variation types
        })
           
        
    
    # Convert results into a DataFrame
    results_df = pd.DataFrame(results)
    
    # Define output path and ensure the directory exists
    output_path = Path('results/gene_variation_summary.csv')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the results to a new CSV file
    results_df.to_csv(output_path, index=False)
    print(f"Results saved to '{output_path}'")
    
    return results_df

# Example usage
if __name__ == "__main__":
    csv_file = 'results/pangenome_summary.csv'  # Replace with the actual file path
    summary = process_gene_variations(csv_file)
    print(summary)

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
        
    Raises:
        FileNotFoundError: If the specified CSV file is not found.
        KeyError: If the expected columns (e.g., 'GENE', 'TYPE', etc.) are missing from the input CSV.
        ValueError: If the CSV file contains unexpected data or formats.
    """
    # Read the CSV file into a pandas DataFrame
    try:
        data = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: The file '{csv_file}' was not found.")
        raise
    except Exception as e:
        print(f"Error: An unexpected error occurred while reading the file '{csv_file}': {e}")
        raise
    
    # Initialize the results list
    results = []
    
    try:
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

            except KeyError as e:
                print(f"Error: Missing expected column '{e.args[0]}' in the input data.")
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
    
    except KeyError as e:
        print(f"Error: Missing expected column '{e.args[0]}' while processing data.")
        raise
    except Exception as e:
        print(f"Error: An unexpected error occurred during processing: {e}")
        raise
    
    # Convert results into a DataFrame
    results_df = pd.DataFrame(results)
    
    # Define output path and ensure the directory exists
    output_path = Path('results/EXTENDED_pangenome_variation_summary.csv')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the results to a new CSV file
    try:
        results_df.to_csv(output_path, index=False)
        print(f"Results saved to '{output_path}'")
    except Exception as e:
        print(f"Error: Failed to save results to '{output_path}': {e}")
        raise
    
    return results_df


def ensure_all_genes_covered(processed_data, reference_csv):
    """
    Ensures all genes from the reference gene list are covered in the processed data.
    If a gene is missing, it will be added with default values. Gene lengths and total lengths
    will be copied from the flank counterpart if available.
    
    Args:
        processed_data (pd.DataFrame): The processed gene variations summary.
        reference_csv (str): Path to the reference gene list.
        
    Returns:
        pd.DataFrame: Updated DataFrame with all genes accounted for.
        
    Raises:
        FileNotFoundError: If the specified reference CSV file is not found.
        KeyError: If the 'GENE' column is missing in the processed data or reference file.
    """
    try:
        reference_data = pd.read_csv(reference_csv)
    except FileNotFoundError:
        print(f"Error: The reference file '{reference_csv}' was not found.")
        raise
    except Exception as e:
        print(f"Error: An unexpected error occurred while reading the reference file '{reference_csv}': {e}")
        raise
    
    try:
        reference_data.rename(columns={'GeneName': 'GENE'}, inplace=True)
        reference_genes = set(reference_data['GENE'])
        
        processed_genes = set(processed_data['GENE'])
        missing_genes = reference_genes - processed_genes

        for gene in missing_genes:
            # Look for the flank counterpart
            flank_gene = f"{gene}_flank"
            flank_data = processed_data.loc[processed_data['GENE'] == flank_gene]

            if not flank_data.empty:
                # Copy lengths from the flank
                gene_length = flank_data.iloc[0]['GENE_LENGTH']
                total_length = flank_data.iloc[0]['TOTAL_LENGTH']
            else:
                gene_length = 'NA'
                total_length = 'NA'
            
            processed_data = pd.concat([
                processed_data,
                pd.DataFrame([{
                    'GENE': gene,
                    'NUM_VARIATIONS': 0,
                    'VARIATION_TYPES': "SNP",  # Default for missing genes
                    'GENE_LENGTH': gene_length,
                    'TOTAL_LENGTH': total_length
                }])
            ])
    
    except KeyError as e:
        print(f"Error: Missing expected column '{e.args[0]}' while ensuring all genes are covered.")
        raise
    except Exception as e:
        print(f"Error: An unexpected error occurred while processing the reference file: {e}")
        raise
    
    # Sort processed data by gene name for consistency
    processed_data = processed_data.sort_values(by='GENE').reset_index(drop=True)
    return processed_data

if __name__ == "__main__":
    csv_file = 'data/EXTENDED_pangenome_summary.csv'  # Replace with actual variation file path
    reference_csv = 'data/SNP-densities-and-RNA.csv'  # Replace with actual reference file path

    try:
        # Step 1: Process variations
        processed_summary = process_gene_variations(csv_file)
        
        # Step 2: Ensure all genes are accounted for
        final_summary = ensure_all_genes_covered(processed_summary, reference_csv)

        # Save final summary to file
        final_output_path = Path('results/EXTENDED_pangenome_variation_summary.csv')
        final_summary.to_csv(final_output_path, index=False)
        print(f"Final results saved to '{final_output_path}'")

    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")

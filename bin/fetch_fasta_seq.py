import os
import csv
import requests
import time

ENSEMBL_REST_API_URL = "https://rest.ensembl.org"

def get_ensembl_id_from_symbol(symbol, species="human"):
    """
    Fetches the Ensembl ID for a given gene symbol using Ensembl REST API.

    Args:
        symbol (str): The gene symbol.
        species (str): The species (default is 'human').

    Returns:
        str: Ensembl ID for the gene.
    """
    url = f"{ENSEMBL_REST_API_URL}/lookup/symbol/{species}/{symbol}?content-type=application/json"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        # Extract Ensembl ID from the response
        data = response.json()
        ensembl_id = data.get('id', None)
        if ensembl_id:
            return ensembl_id
        else:
            raise ValueError(f"Ensembl ID not found for symbol: {symbol}")
    except requests.exceptions.RequestException as e:
        raise ValueError(f"Error fetching Ensembl ID for symbol {symbol}: {e}")

def fetch_fasta_sequence_by_ensembl_id(gene_symbol, ensembl_id, output_dir):
    """
    Fetches the FASTA sequence for a given Ensembl ID from Ensembl REST API.

    Args:
        ensembl_id (str): The Ensembl gene ID.
        output_dir (str): The directory to save the FASTA file.

    Returns:
        str: Path to the saved FASTA file.
    """
    # URL for Ensembl REST API to fetch the FASTA sequence using Ensembl ID
    url = f"{ENSEMBL_REST_API_URL}/sequence/id/{ensembl_id}?content-type=text/x-fasta"
    
    try:
        # Send the GET request to fetch the FASTA sequence for the Ensembl ID
        response = requests.get(url, timeout=10)
        response.raise_for_status()  # Raise an exception for HTTP errors
        
        # If the response contains the FASTA sequence, save it to a file
        fasta_sequence = response.text
        
        # Save the FASTA sequence to a file
        fasta_file_path = os.path.join(output_dir, f"{gene_symbol}.fasta")
        with open(fasta_file_path, "w") as fasta_file:
            fasta_file.write(fasta_sequence)
        
        return fasta_file_path
    except requests.exceptions.RequestException as e:
        raise ValueError(f"Error fetching FASTA sequence for Ensembl ID {ensembl_id}: {e}")

def process_gene_list(csv_file_path, output_dir):
    """
    Reads a CSV file containing gene symbols and locations, then downloads their FASTA sequences.

    Args:
        csv_file_path (str): Path to the CSV file.
        output_dir (str): Directory to save FASTA files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with open(csv_file_path, "r") as csv_file:
        reader = csv.DictReader(csv_file)  # Use DictReader to access columns by name
        for row in reader:
            gene_symbol = row["GeneName"].strip()  # Strip whitespace
            try:
                print(f"Fetching Ensembl ID for Gene Symbol: {gene_symbol}")
                ensembl_id = get_ensembl_id_from_symbol(gene_symbol)  # Get Ensembl ID
                print(f"Fetching FASTA for Ensembl ID: {ensembl_id}")
                fasta_path = fetch_fasta_sequence_by_ensembl_id(gene_symbol, ensembl_id, output_dir)
                print(f"Saved to {fasta_path}")
            except Exception as e:
                print(f"Failed to fetch FASTA for Gene Symbol: {gene_symbol}. Error: {e}")
            
            # Sleep briefly to respect Ensembl's API rate limit (10 requests per second)
            time.sleep(0.1)  # 0.1 seconds = 10 requests per second

if __name__ == "__main__":
    # Input CSV file
    input_csv = "data/SNP-densities-and-RNA.csv"  # Replace with the actual file path

    # Output directory for FASTA files
    output_directory = "data"

    process_gene_list(input_csv, output_directory)

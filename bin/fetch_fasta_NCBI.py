import requests
import os
import csv

# Define the base URL for NCBI E-utilities
base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Chromosome accession mapping (RefSeq accession IDs)
chromosome_accessions = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10"
}

# Function to fetch FASTA sequence using chromosomal range
def fetch_fasta_from_range(chromosome, start, stop):
    """
    Fetches a FASTA sequence from NCBI based on a given chromosomal range.

    Args:
    - chromosome (str): Chromosome identifier (e.g., '1', 'X').
    - start (int): Start position of the range.
    - stop (int): End position of the range.

    Returns:
    - str: FASTA sequence as a string if successful, None otherwise.
    """
    # Check if the chromosome is valid
    if chromosome not in chromosome_accessions:
        print(f"Invalid chromosome: {chromosome}. Skipping.")
        return None

    # Get the RefSeq accession ID for the chromosome
    chromosome_id = chromosome_accessions[chromosome]

    # Construct the URL with chromosomal range parameters
    url = f"{base_url}?db=nuccore&id={chromosome_id}&seq_start={start}&seq_stop={stop}&rettype=fasta"
    
    # Print the URL for debugging purposes
    print(f"Requesting URL: {url}")
    
    try:
        # Make the HTTP request
        response = requests.get(url)
        response.raise_for_status()  # Will raise an HTTPError for bad responses (4xx/5xx)
        
        # Check if the response is successful (status code 200)
        if response.status_code == 200:
            print(f"Successfully fetched data for Chromosome {chromosome}, {start}-{stop}")
            return response.text
        else:
            print(f"Error fetching sequence for Chromosome {chromosome}, {start}-{stop}: {response.status_code}")
            return None

    except requests.exceptions.RequestException as e:
        # Handle any request exceptions (e.g., connection errors)
        print(f"Request error occurred while fetching data for Chromosome {chromosome}, {start}-{stop}: {e}")
        return None

# Function to save the FASTA sequence to a file
def save_fasta(sequence, gene_name):
    """
    Saves the FASTA sequence to a file.

    Args:
    - sequence (str): The FASTA sequence to be saved.
    - gene_name (str): The gene name used for the filename.
    """
    try:
        # Create the data directory if it doesn't exist
        os.makedirs('data', exist_ok=True)

        # Save the FASTA sequence to a file named after the gene
        with open(f"data/{gene_name}.fasta", "w") as file:
            file.write(sequence)
        print(f"FASTA sequence for {gene_name} saved.")

    except OSError as e:
        # Handle file system errors
        print(f"Error saving FASTA sequence for {gene_name}: {e}")

# Function to read gene locations from a CSV file and fetch FASTA sequences
def fetch_fasta_for_genes_from_csv(csv_file):
    """
    Reads gene locations from a CSV file and fetches FASTA sequences for each gene.

    Args:
    - csv_file (str): Path to the CSV file containing gene information (GeneName, Chromosome, Start, End).
    """
    try:
        # Open and read the CSV file
        with open(csv_file, mode='r') as file:
            reader = csv.DictReader(file)

            # Loop through each row in the CSV file
            for row in reader:
                gene_name = row['GeneName']
                chromosome = row['Chromosome'].replace("chr", "")  # Remove 'chr' prefix if present
                start = int(row['Start'])
                stop = int(row['End'])

                # Print info about the gene being fetched
                print(f"Fetching FASTA for {gene_name} on Chromosome {chromosome} from {start} to {stop}...")

                # Fetch the FASTA sequence
                fasta_sequence = fetch_fasta_from_range(chromosome, start, stop)

                # If the sequence is fetched successfully, save it
                if fasta_sequence:
                    save_fasta(fasta_sequence, gene_name)
                else:
                    print(f"Failed to fetch FASTA for {gene_name}.")
    
    except FileNotFoundError as e:
        # Handle file not found error
        print(f"Error: The file '{csv_file}' was not found. Please check the file path.")
    
    except KeyError as e:
        # Handle missing expected columns in CSV
        print(f"Error: Missing expected column in CSV: {e}")
    
    except Exception as e:
        # Catch any other errors
        print(f"An unexpected error occurred: {e}")

# Call the function with your CSV file path
fetch_fasta_for_genes_from_csv("data/SNP-densities-and-RNA.csv")  # Replace with the actual path to your CSV file

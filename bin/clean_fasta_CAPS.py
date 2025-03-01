import re
import os

def process_fasta_file(file_path):
    """
    Processes a FASTA file by:
    1. Converting lowercase nucleotides to uppercase.
    2. Replacing 'Chr' with 'chr'.
    3. Splitting the content into individual chromosomes.
    4. Sorting chromosomes by their numeric identifiers.
    5. Overwriting the original file with the processed content.

    Args:
        file_path (str): Path to the FASTA file to process.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        PermissionError: If there are permission issues with reading or writing the file.
        Exception: For any other unexpected errors during processing.
    """
    try:
        print(f"Starting to process the file: {file_path}")
        
        # Check if file exists
        if not os.path.exists(file_path):
            print(f"Error: The file '{file_path}' does not exist.")
            raise FileNotFoundError(f"File not found: {file_path}")
        
        # Open the file to read its content
        print("Opening file and reading content...")
        with open(file_path, 'r') as file:
            file_content = file.read()
        print(f"File read successfully, length of content: {len(file_content)} characters.")
        
        # Replace lowercase nucleotides with their uppercase counterparts
        print("Processing content: converting lowercase nucleotides to uppercase...")
        processed_content = file_content.translate(str.maketrans('aucgtn', 'AUCGTN'))
        
        # Change 'Chr' to 'chr'
        print("Processing content: changing 'Chr' to 'chr'...")
        processed_content = processed_content.replace('Chr', 'chr')
        
        # Split the content into individual chromosomes
        print("Splitting the content into individual chromosomes...")
        chromosomes = re.split(r'(?=^>chr)', processed_content, flags=re.MULTILINE)
        
        # Sort chromosomes by the numeric part of the chromosome name
        print("Sorting chromosomes by number...")
        sorted_chromosomes = sorted(chromosomes, key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf'))
        
        # Recombine the sorted chromosomes back into the content
        sorted_content = ''.join(sorted_chromosomes)
        
        # Check if any changes were made (this is a simple diagnostic check)
        if processed_content != sorted_content:
            print("Changes detected. Proceeding to overwrite the file.")
        else:
            print("No changes detected, nothing to overwrite.")
        
        # Overwrite the original file with the processed content
        with open(file_path, 'w') as file:
            file.write(sorted_content)
        print(f"File '{file_path}' has been overwritten with processed content.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
        raise
    except PermissionError as e:
        print(f"Error: Permission denied when accessing the file '{file_path}'.")
        raise
    except Exception as e:
        print(f"An unexpected error occurred while processing the file '{file_path}': {e}")
        raise

def process_all_fasta_files(folder_path):
    """
    Processes all FASTA files in a specified folder that start with 'chr' and end with '.fa'.
    
    Args:
        folder_path (str): Path to the folder containing the FASTA files.
    
    Raises:
        ValueError: If no matching files are found in the folder.
    """
    # Get all files in the folder with names starting with 'chr' and ending with '.fa'
    file_list = [f for f in os.listdir(folder_path) if f.startswith('chr') and f.endswith('.fa')]
    
    if not file_list:
        raise ValueError(f"No files found in the folder '{folder_path}' that match the pattern 'chr*.fa'.")
    
    # Process each file
    for file_name in file_list:
        file_path = os.path.join(folder_path, file_name)
        process_fasta_file(file_path)

# Example usage:
folder_path = 'data/Chromosomes_FASTA/'  # Path to the folder containing the chromosome FASTA files

try:
    # process_all_fasta_files(folder_path)
    process_fasta_file("data/hg38.fa")
except Exception as e:
    print(f"Error: The process failed due to: {e}")

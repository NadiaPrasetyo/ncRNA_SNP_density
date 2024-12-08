import os

def combine_sequences(folder_path, output_file):
    """
    Combine sequences from multiple .fa files in a given folder into a single output file.
    The output file will contain the sequences from all .fa files, with empty lines removed.
    
    Args:
        folder_path (str): Path to the folder containing the .fa files.
        output_file (str): Path to the output file where all sequences will be written.
        
    Raises:
        FileNotFoundError: If the folder specified does not exist.
        PermissionError: If there is a permission issue with reading or writing the files.
        ValueError: If no .fa files are found in the specified folder.
    """
    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    
    # Open the output file in write mode
    try:
        with open(output_file, 'w') as outfile:
            # Flag to check if any .fa files are processed
            files_processed = False
            
            # Loop through each file in the folder
            for filename in os.listdir(folder_path):
                if filename.endswith('.fa'):  # Check if the file is a .fa file
                    files_processed = True
                    file_path = os.path.join(folder_path, filename)
                    
                    # Open the .fa file and process its contents
                    try:
                        with open(file_path, 'r') as infile:
                            for line in infile:
                                # Strip leading/trailing whitespaces and skip empty lines
                                stripped_line = line.strip()
                                if stripped_line:  # Only write non-empty lines
                                    outfile.write(stripped_line + '\n')
                    except PermissionError:
                        print(f"Error: Permission denied when trying to read the file '{file_path}'.")
                        raise
                    except Exception as e:
                        print(f"Error: An unexpected error occurred while processing the file '{file_path}': {e}")
                        raise
            
            if not files_processed:
                print(f"Error: No .fa files found in the folder '{folder_path}'.")
                raise ValueError(f"No .fa files found in folder: {folder_path}")
        
        print(f"All sequences have been combined into '{output_file}', with empty lines removed.")
    
    except PermissionError:
        print(f"Error: Permission denied when trying to write to the output file '{output_file}'.")
        raise
    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")
        raise

# Example usage
folder_path = 'data/FASTA'  # Path to the folder containing the .fa files
output_file = 'data/combined_sequences.fa'  # Output file where all sequences will be written

try:
    combine_sequences(folder_path, output_file)
except Exception as e:
    print(f"Error: The process failed due to: {e}")

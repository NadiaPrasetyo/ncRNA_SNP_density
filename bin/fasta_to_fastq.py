import os

def fasta_to_fastq(fasta_dir, fastq_dir):
    """
    Converts all FASTA files in the given directory to FASTQ format.

    Args:
    - fasta_dir (str): The path to the directory containing FASTA files.
    - fastq_dir (str): The path to the directory where FASTQ files will be saved.

    Raises:
    - FileNotFoundError: If the input FASTA directory does not exist.
    - OSError: If there is an issue with creating output FASTQ directory or writing files.
    """
    try:
        # Ensure the output directory exists
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)
            print(f"Created directory: {fastq_dir}")

        # Check if the fasta directory exists
        if not os.path.exists(fasta_dir):
            raise FileNotFoundError(f"The input FASTA directory '{fasta_dir}' does not exist.")

        # Loop through all files in the fasta directory
        for file_name in os.listdir(fasta_dir):
            if file_name.endswith(".fa") or file_name.endswith(".fasta"):
                fasta_file_path = os.path.join(fasta_dir, file_name)
                fastq_file_path = os.path.join(fastq_dir, file_name.replace(".fa", ".fq").replace(".fasta", ".fq"))
                
                try:
                    with open(fasta_file_path, 'r') as fasta_file, open(fastq_file_path, 'w') as fastq_file:
                        header = None
                        sequence = None
                        
                        for line in fasta_file:
                            line = line.strip()

                            # If it's a header line (starts with '>')
                            if line.startswith(">"):
                                if header is not None and sequence is not None:
                                    # Write previous sequence to FASTQ
                                    fastq_file.write(f"@{header}\n")
                                    fastq_file.write(f"{sequence}\n")
                                    fastq_file.write("+\n")
                                    fastq_file.write("I" * len(sequence) + "\n")  # Assign 'I' for quality scores
                                
                                # Update header to the new one and reset sequence
                                header = line[1:]  # Remove the '>' symbol
                                sequence = ""
                            else:
                                # Append to the sequence line
                                sequence += line
                        
                        # Write the last sequence after exiting the loop
                        if header is not None and sequence is not None:
                            fastq_file.write(f"@{header}\n")
                            fastq_file.write(f"{sequence}\n")
                            fastq_file.write("+\n")
                            fastq_file.write("I" * len(sequence) + "\n")  # Assign 'I' for quality scores

                except OSError as e:
                    # Handle file I/O errors (e.g., issues reading or writing files)
                    print(f"Error processing file '{fasta_file_path}': {e}")
                except Exception as e:
                    # Catch any other unexpected errors
                    print(f"An unexpected error occurred while processing '{fasta_file_path}': {e}")

        print(f"FASTA files have been successfully converted to FASTQ format in '{fastq_dir}'.")

    except FileNotFoundError as e:
        print(e)  # Print error if the input directory doesn't exist
    except OSError as e:
        print(f"Error with file system: {e}")  # Handle general file system errors
    except Exception as e:
        print(f"An unexpected error occurred: {e}")  # Catch any unexpected errors


# Example usage:
fasta_to_fastq('data/FASTA', 'data/fastq')

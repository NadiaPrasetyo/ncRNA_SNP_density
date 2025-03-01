import os

def collect_successfully_aligned_reads_with_filenames(folder_path, output_file="data/aligned_reads.sam"):
    """
    Collect all successfully aligned reads (FLAG != 4) from .sam files in a given folder.
    The collected reads are written to an output file, along with the filename they were found in.
    
    Args:
        folder_path (str): Path to the folder containing .sam files.
        output_file (str): Path to the output file for aligned reads with filenames.
        
    Raises:
        FileNotFoundError: If the folder containing .sam files does not exist.
        ValueError: If no .sam files are found in the folder.
        PermissionError: If there are permission issues with reading or writing the files.
        Exception: If any other unexpected error occurs during processing.
    """
    aligned_reads = []
    
    # Ensure the folder exists
    if not os.path.exists(folder_path):
        print(f"Error: The folder '{folder_path}' does not exist.")
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    
    # Check if any .sam files are present
    sam_files_found = False
    
    # Iterate through all files in the directory
    try:
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".sam"):
                sam_files_found = True
                file_path = os.path.join(folder_path, file_name)
                print(f"Processing file: {file_path}")
                
                # Open the SAM file
                try:
                    with open(file_path, "r") as sam_file:
                        for line in sam_file:
                            # Skip header lines
                            if line.startswith("@"):
                                continue
                            
                            # Split the line into fields and check the FLAG
                            fields = line.strip().split("\t")
                            if len(fields) > 1:  # Ensure the line has at least two fields
                                flag = int(fields[1])
                                if flag != 4:  # FLAG != 4 means the read is successfully aligned
                                    # Append the filename and the aligned read to the list
                                    aligned_reads.append(f"{file_name}\t{line}")
                except PermissionError:
                    print(f"Error: Permission denied when trying to read the file '{file_path}'.")
                    raise
                except Exception as e:
                    print(f"Error: An unexpected error occurred while processing the file '{file_path}': {e}")
                    raise
        
        if not sam_files_found:
            print(f"Error: No .sam files found in the folder '{folder_path}'.")
            raise ValueError(f"No .sam files found in folder: {folder_path}")
        
        # Sort the aligned reads alphabetically by the filename (first part before the tab)
        aligned_reads.sort()

        # Write the collected aligned reads with filenames to an output file
        try:
            with open(output_file, "w") as output:
                output.writelines(aligned_reads)
            print(f"Successfully aligned reads with filenames have been saved to {output_file}")
        
        except PermissionError:
            print(f"Error: Permission denied when trying to write to the output file '{output_file}'.")
            raise
        except Exception as e:
            print(f"Error: An unexpected error occurred while writing to the file '{output_file}': {e}")
            raise

    except Exception as e:
        print(f"Error: An unexpected error occurred: {e}")
        raise

# Usage example:
folder_path = "data/Alignments"  # Replace with actual folder path
output_file = "data/aligned_genome_reads.sam"  # Replace with desired output file path

try:
    collect_successfully_aligned_reads_with_filenames(folder_path, output_file)
except Exception as e:
    print(f"Process failed due to error: {e}")

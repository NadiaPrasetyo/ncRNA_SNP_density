import os

def collect_successfully_aligned_reads_with_filenames(folder_path, output_file="data/aligned_reads.sam"):
    """
    Reads all .sam files in a folder, collects lines with FLAG != 4 (successfully aligned reads),
    and records the filename where each read was found.
    
    :param folder_path: Path to the folder containing .sam files
    :param output_file: Path to the output file for aligned reads with filenames
    """
    aligned_reads = []
    
    # Iterate through all files in the directory
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".sam"):
            file_path = os.path.join(folder_path, file_name)
            print(f"Processing file: {file_path}")
            
            # Open the SAM file
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
    
    # Sort the aligned reads alphabetically by the filename (first part before the tab)
    aligned_reads.sort()

    # Write the collected aligned reads with filenames to an output file
    with open(output_file, "w") as output:
        output.writelines(aligned_reads)
    
    print(f"Successfully aligned reads with filenames have been saved to {output_file}")

# Usage example:
# Provide the folder path containing .sam files
folder_path = "data/Alignments"
collect_successfully_aligned_reads_with_filenames(folder_path)

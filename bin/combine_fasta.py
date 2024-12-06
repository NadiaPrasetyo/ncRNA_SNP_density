import os

# Path to the folder containing the .fa files
folder_path = 'data/FASTA'

# Output file where all sequences will be written
output_file = 'data/combined_sequences.fa'

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Loop through each file in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.fa'):  # Check if the file is a .fa file
            file_path = os.path.join(folder_path, filename)
            
            # Open the .fa file and process its contents
            with open(file_path, 'r') as infile:
                for line in infile:
                    # Strip leading/trailing whitespaces and skip empty lines
                    stripped_line = line.strip()
                    if stripped_line:  # Only write non-empty lines
                        outfile.write(stripped_line + '\n')
                
print(f"All sequences have been combined into {output_file}, with empty lines removed.")

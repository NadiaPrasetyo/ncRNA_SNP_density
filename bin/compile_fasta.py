import os

# Define the directory containing the FASTA files
fasta_dir = "data/FASTA"
# Output file where all sequences will be concatenated
output_file = "data/combined_ncRNAs.fasta"

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Iterate through all files in the FASTA directory
    for filename in os.listdir(fasta_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            # Open each FASTA file
            with open(os.path.join(fasta_dir, filename), 'r') as infile:
                # Read and write the content of the current file to the output file
                outfile.write(infile.read())

print(f"All FASTA files have been combined into {output_file}")

import os

# Define the input and output directories
input_dir = 'data/Chromosomes_FASTA'
output_file = 'data/Chromosomes_FASTA/hg38_chromosomes.fa'

# Function to combine FASTA files
def combine_fasta_files(input_directory, output_filename):
    try:
        # Open the output file in write mode
        with open(output_filename, 'w') as outfile:
            # Iterate over all files in the input directory
            for filename in os.listdir(input_directory):
                # Get the full file path
                file_path = os.path.join(input_directory, filename)
                
                # Only process files that have a '.fasta' extension
                if os.path.isfile(file_path) and filename.endswith('.fa'):
                    with open(file_path, 'r') as infile:
                        # Write the content of the file to the output
                        outfile.write(infile.read())
                        # Add a newline between files for clarity
                        outfile.write('\n')
        
        print(f"All FASTA files have been combined into {output_filename}.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Call the function
combine_fasta_files(input_dir, output_file)

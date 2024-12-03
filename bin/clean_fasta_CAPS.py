# Define the path to the input file
file_path = 'data/Chromosomes_FASTA/hg38_chromosomes.fa'

# Function to process the file
def process_fasta_file(file_path):
    print(f"Starting to process the file: {file_path}")
    
    # Open the file to read its content
    print("Opening file and reading content...")
    with open(file_path, 'r') as file:
        file_content = file.read()
    print(f"File read successfully, length of content: {len(file_content)} characters.")
    
    # Replace lowercase nucleotides with their uppercase counterparts
    print("Processing content: converting lowercase nucleotides to uppercase...")
    processed_content = file_content.translate(str.maketrans('aucgt', 'AUCGT'))
    
    # Check if any changes were made (this is a simple diagnostic check)
    if file_content != processed_content:
        print("Changes detected. Proceeding to overwrite the file.")
    else:
        print("No changes detected, nothing to overwrite.")
    
    # Overwrite the original file with the processed content
    with open(file_path, 'w') as file:
        file.write(processed_content)
    print(f"File '{file_path}' has been overwritten with processed content.")

# Run the function to process the file
process_fasta_file(file_path)

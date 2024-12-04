import re

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
    processed_content = file_content.translate(str.maketrans('aucgtn', 'AUCGTN'))
    
    # Change 'Chr' to 'chr'
    print("Processing content: changing 'Chr' to 'chr'...")
    processed_content = processed_content.replace('Chr', 'chr')
    
    # Split the content into individual chromosomes
    print("Splitting the content into individual chromosomes...")
    chromosomes = re.split(r'(?=^>chr)', processed_content, flags=re.MULTILINE)
    
    # Remove duplicates by keeping only the first occurrence of each chromosome
    seen_chromosomes = set()
    unique_chromosomes = []
    
    print("Removing duplicates...")
    for chromosome in chromosomes:
        # Extract chromosome name (e.g., chr1, chr2, etc.)
        match = re.search(r'^>chr\d+', chromosome)
        if match:
            chrom_name = match.group()
            if chrom_name not in seen_chromosomes:
                unique_chromosomes.append(chromosome)
                seen_chromosomes.add(chrom_name)
    
    # Sort chromosomes by the numeric part of the chromosome name
    print("Sorting chromosomes by number...")
    unique_chromosomes = sorted(unique_chromosomes, key=lambda x: int(re.search(r'\d+', x).group()) if re.search(r'\d+', x) else float('inf'))
    
    # Recombine the sorted chromosomes back into the content
    sorted_content = ''.join(unique_chromosomes)
    
    # Check if any changes were made (this is a simple diagnostic check)
    if processed_content != sorted_content:
        print("Changes detected. Proceeding to overwrite the file.")
    else:
        print("No changes detected, nothing to overwrite.")
    
    # Overwrite the original file with the processed content
    with open(file_path, 'w') as file:
        file.write(sorted_content)
    print(f"File '{file_path}' has been overwritten with processed content.")

# Run the function to process the file
process_fasta_file(file_path)

import os

def fasta_to_fastq(fasta_dir, fastq_dir):
    # Ensure the output directory exists
    if not os.path.exists(fastq_dir):
        os.makedirs(fastq_dir)
    
    # Loop through all files in the fasta directory
    for file_name in os.listdir(fasta_dir):
        if file_name.endswith(".fa") or file_name.endswith(".fasta"):
            fasta_file_path = os.path.join(fasta_dir, file_name)
            fastq_file_path = os.path.join(fastq_dir, file_name.replace(".fa", ".fq").replace(".fasta", ".fq"))
            
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

    print(f"FASTA files have been successfully converted to FASTQ format in '{fastq_dir}'.")

# Example usage:
fasta_to_fastq('data/FASTA', 'data/fastq')

import os

def read_bed_files(directory):
    # List all files in the specified directory
    files = [f for f in os.listdir(directory) if f.endswith('.bed')]

    # Process each .bed file
    for file in files:
        file_path = os.path.join(directory, file)
        print(f"\nContents of {file}:")

        try:
            with open(file_path, 'r') as f:
                # Read the header (if any) and first 50 lines
                for i, line in enumerate(f):
                    if i >= 50:
                        break
                    print(line.strip())
        except Exception as e:
            print(f"Error reading {file}: {e}")

# Specify the folder containing the .bed files
folder_path = "data/"

# Call the function to read and print file contents
read_bed_files(folder_path)

#order: chromosome, start, end, name of item, score, strand, coverage/number of reads, percentage of reads that show methylation at this position in the genome
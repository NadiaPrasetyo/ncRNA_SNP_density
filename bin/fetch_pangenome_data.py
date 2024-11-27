
# Define the file path the csv file of genes and locations with high SNP density
csv_file_path = "data/SNP densities and RNA - Sheet1.csv"

# Define the file path to your VCF file
vcf_file_path = "data/decomposed.vcf"

# Specify the number of lines to skip and the number of lines to print
lines_to_skip = 2595
lines_to_print = 10

# Open the VCF file and process it
try:
    with open(vcf_file_path, 'r') as vcf_file:
        # Skip the first 2595 lines
        for _ in range(lines_to_skip):
            if not vcf_file.readline():  # Stop if the end of file is reached early
                print("Reached the end of the file while skipping lines.")
                break
        
        # Print the next 10 lines
        for _ in range(lines_to_print):
            line = vcf_file.readline()
            if not line:  # Stop if the end of file is reached early
                break
            print(line.strip())  # Print each line without trailing newlines
except FileNotFoundError:
    print(f"The file '{vcf_file_path}' does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")
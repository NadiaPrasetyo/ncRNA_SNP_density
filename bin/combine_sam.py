import os

def merge_sam_files(input_folder, output_file):
    # Ensure the input folder exists
    if not os.path.exists(input_folder):
        print(f"Error: The folder '{input_folder}' does not exist.")
        return

    # List all SAM files in the folder
    sam_files = [f for f in os.listdir(input_folder) if f.endswith('.sam')]
    if not sam_files:
        print(f"No SAM files found in the folder '{input_folder}'.")
        return

    print(f"Found {len(sam_files)} SAM files. Merging...")

    headers = set()
    output_lines = []

    for sam_file in sam_files:
        file_path = os.path.join(input_folder, sam_file)
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    # Add headers to the set to ensure uniqueness
                    headers.add(line)
                else:
                    # Add the alignment data to the output list
                    output_lines.append(line)

    # Write the merged content to the output file
    with open(output_file, 'w') as out_f:
        # Write unique headers first
        out_f.writelines(sorted(headers))
        # Write all alignment data
        out_f.writelines(output_lines)

    print(f"Merged SAM file created at '{output_file}'.")

# Example usage
input_folder = "data/Alignments"  # Replace with your folder name
output_file = "data/merged_alignments.sam"  # Replace with your desired output file name
merge_sam_files(input_folder, output_file)

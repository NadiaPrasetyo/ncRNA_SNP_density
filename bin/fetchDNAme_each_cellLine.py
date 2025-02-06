import h5py
import os

# Function to load bigWig file and print first 5 lines
def print_first_lines(file_path):
    try:
        with h5py.File(file_path, 'r') as f:
            print(f"Contents of {file_path}:")
            for i, (key, value) in enumerate(f.items()):
                if i == 5:  # We only want to print the first 5 items
                    break
                print(f"{key}: {value}")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

# List your files
file_directory = 'data/'  # assuming your files are in a folder called 'data'
files = [f for f in os.listdir(file_directory) if f.endswith('.bigWig')]  # Get only BigWig files

# Loop through the files and print the first 5 lines
for file in files:
    print_first_lines(os.path.join(file_directory, file))

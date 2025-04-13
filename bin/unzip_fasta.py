import gzip
import shutil
import os

def unzip_gz_files_recursively(root_directory):
    # Walk through all subdirectories and files
    for dirpath, dirnames, filenames in os.walk(root_directory):
        for filename in filenames:
            if filename.endswith('.gz'):
                gz_path = os.path.join(dirpath, filename)
                output_path = os.path.join(dirpath, filename[:-3])  # Remove .gz

                try:
                    with gzip.open(gz_path, 'rb') as f_in:
                        with open(output_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    print(f'Unzipped: {gz_path} -> {output_path}')
                except Exception as e:
                    print(f'Failed to unzip {gz_path}: {e}')

# Example usage
if __name__ == "__main__":
    target_directory = "data/1000genomes/test/data/NA12282/sequence_read/"  # Replace this with your actual path
    unzip_gz_files_recursively(target_directory)


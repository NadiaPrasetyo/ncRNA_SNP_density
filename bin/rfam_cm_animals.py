import os
import subprocess

# Define paths and file names
data_dir = "data/animals"
rfam_dir = os.path.join(data_dir, "rfam")
output_dir = os.path.join(data_dir, "output")

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# List of genomes and RFAM CMs
genomes = ["galGal6", "danRer11", "rn7", "mm39", "bosTau9"]
rfam_cms = {
    "U8": "RF00096.cm",
    "u5": "RF00020.cm",
    "tRNA": "RF00005.cm",
    "snarA": "RF02556.cm"
}

# Function to calculate Z value for a genome
def calculate_z_value(genome_path):
    try:
        result = subprocess.run(
            ["esl-seqstat", genome_path],
            capture_output=True, text=True, check=True
        )
        print(f"esl-seqstat output for {genome_path}:\n{result.stdout}")  # Debugging line
        for line in result.stdout.split("\n"):
            if line.lstrip().startswith("Total # residues:"):
                try:
                    # Strip extra whitespace and commas, then convert to integer
                    total_residues = int(line.split(":")[1].strip().replace(",", "").replace(" ", ""))
                    z_value = (total_residues * 2) / 1_000_000
                    print(f"Calculated Z value for {genome_path}: {z_value:.6f}")  # Debugging line
                    return z_value
                except ValueError as e:
                    print(f"Error parsing 'Total # of residues' for {genome_path}: {e}")
    except subprocess.CalledProcessError as e:
        print(f"Error calculating Z value for {genome_path}: {e}")
    return None

# Ensure Rfam CM files are pressed
for rfam_name, rfam_cm in rfam_cms.items():
    rfam_path = os.path.join(rfam_dir, rfam_cm)
    try:
        print(f"Pressing CM file: {rfam_path}")
        subprocess.run(["cmpress", rfam_path], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error pressing CM file {rfam_path}: {e}")

# Iterate over genomes and RFAM CMs
for genome in genomes:
    genome_path = os.path.join(data_dir, f"{genome}.fa")
    z_value = calculate_z_value(genome_path)

    if z_value is None:
        print(f"Skipping genome {genome} due to Z value calculation error.")
        continue

    for rfam_name, rfam_cm in rfam_cms.items():
        rfam_path = os.path.join(rfam_dir, rfam_cm)

        # Define output file names
        cmscan_file = os.path.join(output_dir, f"{genome}-{rfam_name}.cmscan")
        tblout_file = f"{cmscan_file}.tblout"
        deoverlapped_tblout_file = f"{cmscan_file}.deoverlapped.tblout"

        # Check if output files already exist
        if os.path.exists(cmscan_file) and os.path.exists(tblout_file) and os.path.exists(deoverlapped_tblout_file):
            print(f"Output files for {genome} with RFAM {rfam_name} already exist. Skipping.")
            continue

        # Construct the cmscan command
        cmscan_cmd = (
            f"cmscan -Z {z_value:.6f} --cut_ga --rfam --nohmmonly "
            f"--tblout {tblout_file} --fmt 2 {rfam_path} {genome_path} > {cmscan_file}"
        )

        # Run the command and capture any errors
        try:
            print(f"Running: {cmscan_cmd}")
            subprocess.run(cmscan_cmd, shell=True, check=True)

            # After cmscan completes, remove overlapping hits based on 'olp' column
            print(f"Removing lower-scoring overlaps from {tblout_file} and saving to {deoverlapped_tblout_file}")
            subprocess.run(f"grep -v ' = ' {tblout_file} > {deoverlapped_tblout_file}", shell=True, check=True)

        except subprocess.CalledProcessError as e:
            print(f"Error running cmscan for genome {genome} with RFAM {rfam_name}: {e}")

print("Processing complete. Results are in the output directory.")
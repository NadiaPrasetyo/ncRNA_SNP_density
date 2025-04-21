import os
import subprocess

# 🔧 Set the range of lines to process here
start_line = 0  # starting from line 6 (0-based index)
end_line =  5  # up to but not including line 11

index_file = "data/1000genomes/1000genomes.sequence.index"
download_dir = "data/1000genomes"
os.makedirs(download_dir, exist_ok=True)

def is_valid_url(url):
    return url.startswith("ftp://") or url.startswith("http://") or url.startswith("https://")

def normalize_ftp_path(path):
    # Ensure only double slashes after the protocol
    return path.replace("ftp:/", "ftp://").replace("ftp:///","ftp://")

with open(index_file, "r") as f:
    lines = [line.strip() for line in f if not line.startswith("#") and line.strip()]

for i, line in enumerate(lines[start_line:end_line], start=start_line):
    fields = line.split('\t')
    if len(fields) < 20:
        print(f"Skipping malformed line {i + 1}")
        continue

    primary_url = normalize_ftp_path(fields[0])
    paired_url = normalize_ftp_path(fields[19]) if len(fields) > 19 else ""

    for url in [primary_url, paired_url]:
        if not url or not is_valid_url(url):
            print(f"❌ Invalid or missing URL: {url}")
            continue

        filename = os.path.basename(url)
        dest_path = os.path.join(download_dir, filename)

        print(f"Downloading {filename} to {download_dir}...")
        try:
            subprocess.run(["wget", "-O", dest_path, url], check=True)
            print(f"✅ Downloaded: {url}")
        except subprocess.CalledProcessError:
            print(f"❌ Failed to download: {url}")

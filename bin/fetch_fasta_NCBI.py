import requests
import os

# List of gene symbols (replace with your actual gene list)
gene_symbols = [
    "SNAR-B2", "TRG-CCC5-1", "SNAR-C4", "TRK-TTT3-2", "TRE-TTC3-1", "TRG1", 
    "TRA-CGC5-1-001", "tRNA-Gly-CCC-4-1", "tRNA-Glu-CTC-1-2", "tRNA-Ile-AAT-12-1",
    "tRNA-Gly-GCC-2-4", "tRNA-Pro-AGG-2-4", "TRG-CCC6-1", "TRA-TGC4-1", "TRP-AGG1-1", 
    "TRK-CTT3-1", "tRNA-Glu-TTC-4-1", "tRNA-Thr-CGT-5-1", "TRG-TCC1-1", "tRNA-Val-CAC-3-1", 
    "tRNA-Thr-TGT-4-1", "TRV-TAC1-1", "SNAR-A1", "tRNA-Pro-TGG-1-1", "tRNA-Gln-CTG-2-1", 
    "tRNA-Arg-CCT-3-1", "TRV", "TRF-GAA1-4", "tRNA-iMet-CAT-2-1", "SNAR-C3", "TRA-AGC11-1", 
    "TRL-TAG2-1", "tRNA-Met-CAT-6-1", "TRM-CAT3-2", "tRNA-Leu-CAG-2-2", "TRA-AGC1-1", 
    "tRNA-Gly-TCC-3-1", "tRNA-Lys-CTT-5-1", "tRNA-Ile-AAT-2-1", "TRA-CGC4-1", "TRK-TTT4-1", 
    "SNAR-G2", "tRNA-Asn-GTT-1-1", "TRE-TTC2-1", "tRNA-Leu-TAG-3-1", "Val-tRNA", "TRA-AGC20-1", 
    "tRNA-Pro-TGG-3-3", "tRNA-Arg-ACG-1-2", "TRM-CAT2-1", "tRNA-Thr-CGT-6-1", "TRA-AGC5-1", 
    "tRNA-Val-CAC-2-1", "tRNA-Lys-CTT-1-2", "TRI-AAT4-1", "tRNA-Lys-TTT-1-1", "TRA-AGC4-1", 
    "tRNA-Lys-CTT-2-1", "TRI-AAT5-1", "tRNA-Thr-TGT-1-1", "tRNA-Asp-GTC-1-1", "tRNA-Lys-TTT-6-1", 
    "TRF-GAA2-1", "TRP-TGG2-1", "TRT-TGT5-1", "TRR-TCG1-1", "tRNA-Ala-CGC-3-1", "TRE-TTC5-1", 
    "TRV-AAC5-1", "TRR-ACG2-3", "TRT-AGT1-2", "TRL-TAA4-1", "tRNA-Ala-TGC-3-1", "TRR-TCG3-1", 
    "TRS-AGA2-6", "tRNA-Ser-TGA-2-1", "lnc-SLCO4A1-8", "tRNA-Ile-AAT-1-1", "tRNA-Gln-TTG-3-1", 
    "TRQ-CTG6-1", "TRK-TTT7-1", "TRT-TGT2-1", "TRV-CAC1-5", "tRNA-Ala-CGC-1-1", "TRX-CAT1-5", 
    "tRNA-Val-CAC-1-3", "SNAR-G1", "tRNA-Arg-TCG-4-1", "TRC-GCA11-1", "TRD-GTC2-8", "TRA-AGC15-1",
    "tRNA-Ser-GCT-2-1", "tRNA-Ser", "TRC-GCA4-1", "TRL-AAG2-4", "TRN-GTT2-4", "TRA-CGC2-1", 
    "TRG-CCC2-1", "TRQ-TTG2-1", "tRNA-Val-AAC-3-1", "tRNA-Ala-TGC-2-1", "tRNA-Asn-GTT-4-1", 
    "tRNA-Ser-CGA-2-1", "tRNA-Asn-GTT-3-2", "SNAR-F", "tRNA-Lys-TTT-5-1", "tRNA-Arg-CCT-2-1", 
    "TRS-AGA2-6"
]

# NCBI Gene REST API base URL
base_url = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha2/gene"

# Function to fetch FASTA sequence for a gene symbol
def fetch_fasta_by_symbol(gene_symbol):
    # Construct the API URL for downloading the gene data (FASTA format)
    url = f"{base_url}/symbol/{gene_symbol}/download"
    
    response = requests.get(url)
    
    # If the response is successful, return the FASTA data
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error fetching data for {gene_symbol}: {response.status_code}")
        return None

# Function to save the FASTA sequence to a file
def save_fasta(sequence, gene_symbol):
    # Create the data directory if it doesn't exist
    os.makedirs('data', exist_ok=True)
    
    # Save the FASTA sequence to a file named after the gene symbol
    with open(f"data/{gene_symbol}.fasta", "w") as file:
        file.write(sequence)
    print(f"FASTA sequence for {gene_symbol} saved.")

# Loop through the gene symbols and fetch their FASTA sequences
for gene_symbol in gene_symbols:
    print(f"Fetching FASTA for {gene_symbol}...")
    fasta_sequence = fetch_fasta_by_symbol(gene_symbol)
    
    if fasta_sequence:
        save_fasta(fasta_sequence, gene_symbol)
    else:
        print(f"Failed to fetch FASTA for {gene_symbol}.")
    
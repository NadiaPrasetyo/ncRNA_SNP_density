import os
import json
import csv

# Directory containing the JSON files
data_folder = 'data/BLASTN_RefSeq'
path = 'data/BLASTn_results_summary.csv'

# Function to extract the desired information from a single JSON file
def extract_info_from_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Check if the expected top-level key exists
    if 'BlastOutput2' not in data:
        print(f"Warning: 'BlastOutput2' not found in {file_path}. Skipping this file.")
        return None  # Return None if the file doesn't have the expected structure

    # Extracting query_title
    query_title = data['BlastOutput2']['report']['results']['search']['query_title']
    
    # Dictionary to store unique hit titles and their HSP data
    hits_info = {}

    for hit in data['BlastOutput2']['report']['results']['search']['hits']:
        for description in hit['description']:
            hit_title = description['title']  # Extracting the description: title
            
            # Filter: Only include hits with "Primary Assembly" in the hit title
            if "Primary Assembly" in hit_title:
                # Initialize the hit in the dictionary if not already present
                if hit_title not in hits_info:
                    hits_info[hit_title] = {
                        'max_score': -float('inf'),
                        'total_score': 0,
                        'query_from': float('inf'),
                        'query_to': -float('inf'),
                        'e_value': float('inf'),
                        'identity': 0,
                        'align_len': 0,
                        'accession': description.get('accession', 'Unknown')
                    }

                # Extracting hsps (high-scoring pairs)
                for hsp in hit['hsps']:
                    # Update max score and total score
                    hits_info[hit_title]['max_score'] = max(hits_info[hit_title]['max_score'], hsp['bit_score'])
                    hits_info[hit_title]['total_score'] += hsp['score']
                    
                    # Update query coverage
                    hits_info[hit_title]['query_from'] = min(hits_info[hit_title]['query_from'], hsp['query_from'])
                    hits_info[hit_title]['query_to'] = max(hits_info[hit_title]['query_to'], hsp['query_to'])
                    
                    # Update e-value (taking the smallest)
                    hits_info[hit_title]['e_value'] = min(hits_info[hit_title]['e_value'], hsp['evalue'])
                    
                    # Update identity and alignment length for percentage identity
                    hits_info[hit_title]['identity'] += hsp['identity']
                    hits_info[hit_title]['align_len'] += hsp['align_len']

    return query_title, hits_info

# Function to process all JSON files in the directory and extract information
def process_all_files(data_folder):
    results = []
    for filename in os.listdir(data_folder):
        if filename.endswith(".json"):
            file_path = os.path.join(data_folder, filename)
            extracted_data = extract_info_from_json(file_path)
            if extracted_data:  # Only add if data is extracted
                query_title, hits_info = extracted_data
                results.append({
                    'query_title': query_title,
                    'hits': hits_info
                })
    return results

# Function to calculate query cover percentage
def calculate_query_cover(query_from, query_to, query_length):
    return (query_to - query_from + 1) / query_length * 100

# Process the files and prepare data for CSV output
extracted_data = process_all_files(data_folder)

# Open the CSV file to write the results
with open(path, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header row
    writer.writerow(['GENE', 'Hit Title', 'Scientific Name', 'Max Score', 'Total Score', 'Query Cover', 'E value', 'Per. ident', 'Acc. Len'])
    
    if extracted_data:
        for entry in extracted_data:
            query_title = entry['query_title']  # GENE field will be the query title
            for hit_title, hit_data in entry['hits'].items():
                # Calculate query cover and percentage identity
                query_cover = calculate_query_cover(hit_data['query_from'], hit_data['query_to'], len(entry['query_title']))
                percent_identity = (hit_data['identity'] / hit_data['align_len']) * 100 if hit_data['align_len'] else 0

                # Prepare the summary data
                scientific_name = 'Homo sapiens'  # You can adjust this if needed
                max_score = hit_data['max_score']
                total_score = hit_data['total_score']
                e_value = hit_data['e_value']
                align_len = hit_data['align_len']

                # Write row to CSV
                writer.writerow([query_title, hit_title, scientific_name, max_score, total_score, f'{query_cover:.2f}%', e_value, f'{percent_identity:.2f}', align_len])
    else:
        print("No valid data extracted from the files.")

print(f"Results have been written to {path}")

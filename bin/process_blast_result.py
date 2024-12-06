import os
import json
import csv

# Directory containing the JSON files
data_folder = 'data/BLASTN_RefSeq'
path = 'data/BLASTn_results_summary.csv'

# Function to extract the desired information from a single JSON file
def extract_info_from_json(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
    except json.JSONDecodeError:
        print(f"Error: Failed to decode JSON in {file_path}. Skipping this file.")
        return None
    except Exception as e:
        print(f"Error: {e} occurred while processing {file_path}. Skipping this file.")
        return None

    # Check if the expected top-level key exists
    if 'BlastOutput2' not in data:
        print(f"Warning: 'BlastOutput2' not found in {file_path}. Skipping this file.")
        return None  # Return None if the file doesn't have the expected structure

    # Extracting query_title and query_length
    search_data = data['BlastOutput2']['report']['results']['search']
    query_title = search_data['query_title']
    query_len = search_data['query_len']  # Query length is directly provided

    # Dictionary to store unique hit titles and their HSP data
    hits_info = {}

    for hit in search_data['hits']:
        for description in hit['description']:
            hit_title = description['title']  # Extracting the description: title
            
            # Filter: Only include hits with "Primary Assembly" in the hit title
            if "Primary Assembly" in hit_title:
                # Initialize the hit in the dictionary if not already present
                if hit_title not in hits_info:
                    hits_info[hit_title] = {
                        'max_score': -float('inf'),
                        'total_score': 0,
                        'max_identity': 0,
                        'max_align_len': 0,
                        'e_value': float('inf'),
                        'accession': description.get('accession', 'Unknown'),
                        'mapped_position': None  # To store the hit_from for the max identity
                    }

                # Extracting hsps (high-scoring pairs)
                for hsp in hit['hsps']:
                    # Update max score and total score
                    hits_info[hit_title]['max_score'] = max(hits_info[hit_title]['max_score'], hsp['bit_score'])
                    hits_info[hit_title]['total_score'] += hsp['bit_score']  # Total score is based on bit_score
                    
                    # Update maximum alignment length
                    hits_info[hit_title]['max_align_len'] = max(hits_info[hit_title]['max_align_len'], hsp['align_len'])
                    
                    # Update e-value (taking the smallest)
                    hits_info[hit_title]['e_value'] = min(hits_info[hit_title]['e_value'], hsp['evalue'])

                    # Check if this HSP has the max identity so far and update
                    if hsp['identity'] > hits_info[hit_title]['max_identity']:
                        hits_info[hit_title]['max_identity'] = hsp['identity']
                        hits_info[hit_title]['mapped_position'] = hsp['hit_from']

    return query_title, query_len, hits_info

# Function to process all JSON files in the directory and extract information
def process_all_files(data_folder):
    results = []
    for filename in os.listdir(data_folder):
        if filename.endswith(".json"):
            file_path = os.path.join(data_folder, filename)
            extracted_data = extract_info_from_json(file_path)
            if extracted_data:  # Only add if data is extracted
                query_title, query_len, hits_info = extracted_data
                results.append({
                    'query_title': query_title,
                    'query_len': query_len,
                    'hits': hits_info
                })
    return results

# Process the files and prepare data for CSV output
extracted_data = process_all_files(data_folder)

# Open the CSV file to write the results
with open(path, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write header row
    writer.writerow(['GENE', 'Hit Title', 'Scientific Name', 'Max Score', 'Total Score', 'Query Cover', 'E value', 'Per. ident', 'Mapped Position'])
    
    if extracted_data:
        for entry in extracted_data:
            query_title = entry['query_title']  # GENE field will be the query title
            query_len = entry['query_len']  # Query length from the JSON data
            for hit_title, hit_data in entry['hits'].items():
                # Calculate query cover as the maximum alignment length percentage
                query_cover = (hit_data['max_align_len'] / query_len) * 100 if query_len else 0
                
                if query_cover > 100:
                    query_cover = 100
                
                # Calculate percentage identity as the maximum identity percentage
                percent_identity = (hit_data['max_identity'] / query_len) * 100 if query_len else 0

                # Prepare the summary data
                scientific_name = 'Homo sapiens'  # You can adjust this if needed
                max_score = hit_data['max_score']
                total_score = hit_data['total_score']
                e_value = hit_data['e_value']
                mapped_position = hit_data['mapped_position']

                # Write row to CSV
                writer.writerow([query_title, hit_title, scientific_name, max_score, total_score, f'{query_cover:.2f}%', e_value, f'{percent_identity:.2f}', mapped_position])
    else:
        print("No valid data extracted from the files.")

print(f"Results have been written to {path}")

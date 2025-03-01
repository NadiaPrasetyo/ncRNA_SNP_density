import os
import json
import csv

# Directory containing the JSON files
data_folder = 'data/BLASTN_RefSeq'
output_csv_path = 'data/BLASTn_results_summary.csv'

def extract_info_from_json(file_path):
    """
    Extracts information from a JSON file produced by BLASTN.

    Parameters:
    file_path (str): Path to the JSON file.

    Returns:
    tuple: (query_title, query_len, hits_info) if successful, otherwise None.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
    except json.JSONDecodeError:
        print(f"Error: Failed to decode JSON in {file_path}. Skipping this file.")
        return None
    except Exception as e:
        print(f"Error: {e} occurred while processing {file_path}. Skipping this file.")
        return None

    if 'BlastOutput2' not in data:
        print(f"Warning: 'BlastOutput2' not found in {file_path}. Skipping this file.")
        return None

    try:
        search_data = data['BlastOutput2']['report']['results']['search']
        query_title = search_data['query_title']
        query_len = search_data['query_len']
        
        hits_info = {}

        for hit in search_data['hits']:
            for description in hit['description']:
                hit_title = description['title']
                if "Primary Assembly" in hit_title:
                    if hit_title not in hits_info:
                        hits_info[hit_title] = {
                            'max_score': -float('inf'),
                            'total_score': 0,
                            'max_identity': 0,
                            'max_align_len': 0,
                            'e_value': float('inf'),
                            'accession': description.get('accession', 'Unknown'),
                            'mapped_position': None
                        }
                    for hsp in hit['hsps']:
                        hits_info[hit_title]['max_score'] = max(hits_info[hit_title]['max_score'], hsp['bit_score'])
                        hits_info[hit_title]['total_score'] += hsp['bit_score']
                        hits_info[hit_title]['max_align_len'] = max(hits_info[hit_title]['max_align_len'], hsp['align_len'])
                        hits_info[hit_title]['e_value'] = min(hits_info[hit_title]['e_value'], hsp['evalue'])
                        if hsp['identity'] > hits_info[hit_title]['max_identity']:
                            hits_info[hit_title]['max_identity'] = hsp['identity']
                            hits_info[hit_title]['mapped_position'] = hsp['hit_from']
        return query_title, query_len, hits_info

    except KeyError as e:
        print(f"Error: Missing expected key {e} in {file_path}. Skipping this file.")
        return None

def process_all_files(data_folder):
    """
    Processes all JSON files in a folder to extract BLASTN results.

    Parameters:
    data_folder (str): Path to the folder containing JSON files.

    Returns:
    list: A list of dictionaries containing query_title, query_len, and hits_info.
    """
    results = []
    for filename in os.listdir(data_folder):
        if filename.endswith(".json"):
            file_path = os.path.join(data_folder, filename)
            extracted_data = extract_info_from_json(file_path)
            if extracted_data:
                query_title, query_len, hits_info = extracted_data
                results.append({
                    'query_title': query_title,
                    'query_len': query_len,
                    'hits': hits_info
                })
    return results

def write_to_csv(data, output_csv_path):
    """
    Writes extracted BLASTN data to a CSV file.

    Parameters:
    data (list): Extracted data from BLASTN JSON files.
    output_csv_path (str): Path to the output CSV file.
    """
    with open(output_csv_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['GENE', 'Hit Title', 'Scientific Name', 'Max Score', 'Total Score', 'Query Cover', 'E value', 'Per. ident', 'Mapped Position'])

        if data:
            for entry in data:
                query_title = entry['query_title']
                query_len = entry['query_len']
                for hit_title, hit_data in entry['hits'].items():
                    query_cover = (hit_data['max_align_len'] / query_len) * 100 if query_len else 0
                    query_cover = min(query_cover, 100)  # Ensure query_cover does not exceed 100
                    percent_identity = (hit_data['max_identity'] / query_len) * 100 if query_len else 0
                    scientific_name = 'Homo sapiens'  # Default value, adjust as needed
                    writer.writerow([
                        query_title, hit_title, scientific_name,
                        hit_data['max_score'], hit_data['total_score'],
                        f'{query_cover:.2f}%', hit_data['e_value'],
                        f'{percent_identity:.2f}', hit_data['mapped_position']
                    ])
        else:
            print("No valid data extracted from the files.")

if __name__ == "__main__":
    try:
        # Process the JSON files and write to CSV
        extracted_data = process_all_files(data_folder)
        write_to_csv(extracted_data, output_csv_path)
        print(f"Results have been written to {output_csv_path}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

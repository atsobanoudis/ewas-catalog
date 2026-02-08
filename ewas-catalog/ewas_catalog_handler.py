import csv

# File paths
cpg_file_path = 'cpgs.txt'
results_file_path = 'ewas-catalog/ewascatalog-results.txt'
studies_file_path = 'ewas-catalog/ewascatalog-studies.txt'
output_file_path = 'ewas-catalog/combined_ewas_results.csv'

print("Starting combined extraction...")

# 1. Read the list of target CpGs
with open(cpg_file_path, 'r') as f:
    target_cpgs = set(line.strip() for line in f if line.strip())
print(f"Loaded {len(target_cpgs)} target CpGs.")

# 2. Load studies into a dictionary for quick lookup
studies_data = {}
with open(studies_file_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        # Map StudyID to the whole row
        studies_data[row['StudyID']] = row
print(f"Loaded {len(studies_data)} studies.")

# 3. Process results and merge with study details
with open(results_file_path, 'r', encoding='utf-8') as f_in, \
     open(output_file_path, 'w', newline='', encoding='utf-8') as f_out:
    
    reader = csv.DictReader(f_in, delimiter='\t')
    
    # Define combined fieldnames
    # Take study fieldnames from the first study found
    sample_study = next(iter(studies_data.values()))
    study_fieldnames = [f for f in sample_study.keys() if f != 'StudyID']
    fieldnames = reader.fieldnames + study_fieldnames
    
    writer = csv.DictWriter(f_out, fieldnames=fieldnames)
    writer.writeheader()

    match_count = 0
    for row in reader:
        if row['CpG'] in target_cpgs:
            # Get study details for this row's StudyID
            study_info = studies_data.get(row['StudyID'], {})
            
            # Combine the result row with its study info (excluding duplicate StudyID)
            combined_row = {**row}
            for k, v in study_info.items():
                if k != 'StudyID':
                    combined_row[k] = v
            
            writer.writerow(combined_row)
            match_count += 1

print(f"Finished! Saved {match_count} combined records to '{output_file_path}'.")
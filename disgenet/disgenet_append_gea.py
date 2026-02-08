import pandas as pd
import numpy as np

def main():
    # File paths
    gea_path = 'data/disgenet_gea.csv'
    genes_path = 'annotated_genes.xlsx'
    output_path = 'annotated_genes_output.xlsx'

    print("Loading data...")
    try:
        df_gea = pd.read_csv(gea_path)
        df_genes = pd.read_excel(genes_path)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return

    # 1. Subset entire disgenet_gea.csv
    # Criteria: "Mental or Behavioral Dysfunction (T048)" in "diseaseClasses_UMLS_ST"
    target_class = "Mental or Behavioral Dysfunction (T048)"
    
    # Check if column exists
    if 'diseaseClasses_UMLS_ST' not in df_gea.columns:
        print("Error: Column 'diseaseClasses_UMLS_ST' not found in CSV.")
        return

    # Subset
    # We strip whitespace just in case
    df_psych = df_gea[df_gea['diseaseClasses_UMLS_ST'].astype(str).str.strip() == target_class].copy()
    
    print(f"Original GEA rows: {len(df_gea)}")
    print(f"Subsetted Psych rows: {len(df_psych)}")

    if len(df_psych) == 0:
        print("Warning: No rows matched the filter criteria. Output columns will be empty.")

    # Prepare lists for new columns
    psych_diseases_col = []
    evidence_col = []

    # Polarity sorting map
    polarity_order = {
        'Positive': 0,
        'NAPolarity': 1,
        'Negative': 2
    }

    print("Processing genes...")
    
    # Iterate through each gene in annotated_genes.xlsx
    # Assuming 'symbol' is the column name for gene symbol based on previous inspection
    if 'symbol' not in df_genes.columns:
        print("Error: 'symbol' column not found in annotated_genes.xlsx")
        return

    for index, row in df_genes.iterrows():
        gene_symbol = row['symbol']
        
        # Find matching rows in the subsetted GEA data
        # Assuming 'gene_symbol' is the column in GEA data
        gene_data = df_psych[df_psych['gene_symbol'] == gene_symbol]

        if gene_data.empty:
            psych_diseases_col.append(np.nan)
            evidence_col.append(np.nan)
            continue

        # --- Task 2 part 1: Create 'disgenet_psych_diseases' ---
        # "group their 'score' column values. perform a quick check if they all agree"
        
        unique_diseases = gene_data['disease_name'].unique()
        disease_score_map = {} # Store disease -> score (for sorting evidence later too)
        disease_strings = []
        
        valid_diseases = []

        for disease in unique_diseases:
            d_rows = gene_data[gene_data['disease_name'] == disease]
            scores = d_rows['score'].unique()
            
            if len(scores) > 1:
                # User said: "if they don't agree write 'error'"
                # But we need a single score for sorting. 
                # We will mark it as error in the text but maybe use the max or first for sorting?
                # The prompt implies the final output string contains the error or the score.
                # "schizophrenia would all be 0.7... if they don't agree write 'error'"
                final_score_str = "error"
                # For sorting purposes, we'll take the max score so it floats to top if important
                sort_score = d_rows['score'].max() 
            else:
                final_score_str = str(scores[0])
                sort_score = scores[0]
            
            disease_score_map[disease] = sort_score
            valid_diseases.append({
                'name': disease,
                'score_val': sort_score,
                'display': f"{disease}, {final_score_str}"
            })

        # Order by score (descending)
        # "schizophrenia would supercede bipolar disorder" (0.7 vs 0.5)
        valid_diseases.sort(key=lambda x: x['score_val'], reverse=True)
        
        psych_diseases_str = ";\n".join([d['display'] for d in valid_diseases])
        psych_diseases_col.append(psych_diseases_str)

        # --- Task 3: Create 'disgenet_evidence' ---
        # "list each disease name, polarity, pmYear, and reference number"
        # "order disease names together and by score... within each disease subset order them nominally by polarity... further... order them chronologically"
        
        evidence_lines = []
        
        # We can iterate through the sorted diseases from Task 1 to maintain Disease order
        for d_obj in valid_diseases:
            disease_name = d_obj['name']
            d_rows = gene_data[gene_data['disease_name'] == disease_name].copy()
            
            # Helper for polarity sorting
            # Handle NA polarity -> "NAPolarity"
            d_rows['polarity_filled'] = d_rows['polarity'].fillna('NAPolarity')
            
            # Helper for year sorting
            # We treat NaN years as usually last or first? Prompt: "2006, 2007, etc.." implies ascending.
            # We'll put NaN years at the end or keep them as is. 
            # Let's fill NaN with a high number or low number? Usually NA year means unknown.
            # Let's assume 9999 for sorting if NA, or 0.
            # "2006, 2007" -> Ascending.
            d_rows['pmYear_filled'] = d_rows['pmYear'].fillna(9999) 

            # Sort within disease:
            # 1. Polarity (Positive, NAPolarity, Negative)
            # 2. Year (Ascending)
            
            # Create a sorting key for polarity
            d_rows['pol_rank'] = d_rows['polarity_filled'].map(lambda x: polarity_order.get(x, 3)) # 3 for unexpected
            
            d_rows_sorted = d_rows.sort_values(by=['pol_rank', 'pmYear_filled'], ascending=[True, True])
            
            for _, row_ev in d_rows_sorted.iterrows():
                # Format: Disease, Polarity, Year, Ref
                
                # Polarity
                pol = row_ev['polarity']
                if pd.isna(pol):
                    pol_str = "NAPolarity"
                else:
                    pol_str = str(pol)
                
                # Year
                year = row_ev['pmYear']
                if pd.isna(year):
                    year_str = "NA" # Or empty? Prompt says "2004, 1348843". If NA year, maybe just "NA"?
                    # User example: "Schizophrenia, NAPolarity, 2004, 1348843"
                    # User example 2: "Schizophrenia, Positive, 2008, NA(CLINVAR, GeneticVariation)"
                    # Doesn't explicitly show NA year example. I will use "NA".
                else:
                    try:
                        year_str = str(int(year))
                    except ValueError:
                        year_str = str(year)
                
                # Reference
                # "IF reference_type = PMID (if not, then write the reference value (likely will be NA) and include source and assocationType in parenthesis"
                ref_type = str(row_ev['reference_type'])
                ref_val = row_ev['reference']
                source = row_ev['source']
                assoc = row_ev['associationType']
                
                if ref_type == 'PMID':
                    # Sometimes reference is float, convert to int str
                    try:
                        ref_str = str(int(ref_val))
                    except (ValueError, TypeError):
                        ref_str = str(ref_val)
                else:
                    # reference value (likely NA)
                    if pd.isna(ref_val):
                        r_v = "NA"
                    else:
                        r_v = str(ref_val)
                    
                    ref_str = f"{r_v} ({source}, {assoc})"
                
                line = f"{disease_name}, {pol_str}, {year_str}, {ref_str}"
                evidence_lines.append(line)
        
        evidence_col.append(";\n".join(evidence_lines))

    # Add columns to dataframe
    df_genes['disgenet_psych_diseases'] = psych_diseases_col
    df_genes['disgenet_evidence'] = evidence_col

    # Save with formatting
    print(f"Saving to {output_path}...")
    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        df_genes.to_excel(writer, index=False, sheet_name='Sheet1')
        worksheet = writer.sheets['Sheet1']
        
        from openpyxl.styles import Alignment
        wrap_style = Alignment(wrapText=True, vertical='top')
        
        # Find column indices
        cols_to_wrap = ['disgenet_psych_diseases', 'disgenet_evidence']
        for col_name in cols_to_wrap:
            if col_name in df_genes.columns:
                col_idx = df_genes.columns.get_loc(col_name) + 1 # 1-based
                for row in range(2, len(df_genes) + 2): # Start from row 2 (skip header)
                    cell = worksheet.cell(row=row, column=col_idx)
                    cell.alignment = wrap_style

    print("Done.")

if __name__ == "__main__":
    main()

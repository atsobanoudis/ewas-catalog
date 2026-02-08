import pandas as pd
import numpy as np

def annotate_with_disgenet(genes_df: pd.DataFrame, disgenet_df: pd.DataFrame, psych_only: bool = True) -> pd.DataFrame:
    """
    Annotates a gene DataFrame with DisGeNET associations.
    
    Args:
        genes_df: DataFrame containing a 'symbol' column.
        disgenet_df: DataFrame from DisGeNET GEA export.
        psych_only: If True, filters for 'Mental or Behavioral Dysfunction (T048)'.
    """
    df_genes = genes_df.copy()
    
    if psych_only:
        target_class = "Mental or Behavioral Dysfunction (T048)"
        if 'diseaseClasses_UMLS_ST' in disgenet_df.columns:
            df_gea_subset = disgenet_df[disgenet_df['diseaseClasses_UMLS_ST'].astype(str).str.strip() == target_class].copy()
        else:
            print("Warning: diseaseClasses_UMLS_ST column missing, cannot filter for psych traits.")
            df_gea_subset = disgenet_df.copy()
        col_suffix = "_psych"
    else:
        df_gea_subset = disgenet_df.copy()
        col_suffix = ""

    disease_col_name = f"disgenet{col_suffix}_diseases"
    evidence_col_name = f"disgenet{col_suffix}_evidence"

    diseases_col = []
    evidence_col = []

    polarity_order = {
        'Positive': 0,
        'NAPolarity': 1,
        'Negative': 2
    }

    for _, row in df_genes.iterrows():
        gene_symbol = row['symbol']
        gene_data = df_gea_subset[df_gea_subset['gene_symbol'] == gene_symbol]

        if gene_data.empty:
            diseases_col.append(np.nan)
            evidence_col.append(np.nan)
            continue

        unique_diseases = gene_data['disease_name'].unique()
        valid_diseases = []

        for disease in unique_diseases:
            d_rows = gene_data[gene_data['disease_name'] == disease]
            scores = d_rows['score'].unique()
            
            if len(scores) > 1:
                final_score_str = "error"
                sort_score = d_rows['score'].max() 
            else:
                final_score_str = str(scores[0])
                sort_score = scores[0]
            
            valid_diseases.append({
                'name': disease,
                'score_val': sort_score,
                'display': f"{disease}, {final_score_str}"
            })

        valid_diseases.sort(key=lambda x: x['score_val'], reverse=True)
        diseases_str = ";\n".join([d['display'] for d in valid_diseases])
        diseases_col.append(diseases_str)

        evidence_lines = []
        for d_obj in valid_diseases:
            disease_name = d_obj['name']
            d_rows = gene_data[gene_data['disease_name'] == disease_name].copy()
            
            d_rows['polarity_filled'] = d_rows['polarity'].fillna('NAPolarity')
            d_rows['pmYear_filled'] = d_rows['pmYear'].fillna(9999) 
            d_rows['pol_rank'] = d_rows['polarity_filled'].map(lambda x: polarity_order.get(x, 3))
            
            d_rows_sorted = d_rows.sort_values(by=['pol_rank', 'pmYear_filled'], ascending=[True, True])
            
            for _, row_ev in d_rows_sorted.iterrows():
                pol = str(row_ev['polarity']) if pd.notna(row_ev['polarity']) else "NAPolarity"
                
                year = row_ev['pmYear']
                year_str = str(int(year)) if pd.notna(year) else "NA"
                
                ref_type = str(row_ev['reference_type'])
                ref_val = row_ev['reference']
                if ref_type == 'PMID':
                    try:
                        ref_str = str(int(ref_val))
                    except (ValueError, TypeError):
                        ref_str = str(ref_val)
                else:
                    r_v = str(ref_val) if pd.notna(ref_val) else "NA"
                    ref_str = f"{r_v} ({row_ev['source']}, {row_ev['associationType']})"
                
                line = f"{disease_name}, {pol}, {year_str}, {ref_str}"
                evidence_lines.append(line)
        
        evidence_col.append(";\n".join(evidence_lines))

    df_genes[disease_col_name] = diseases_col
    df_genes[evidence_col_name] = evidence_col
    
    return df_genes
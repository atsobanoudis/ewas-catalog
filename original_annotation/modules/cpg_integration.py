import pandas as pd
import numpy as np
import re
import json

def expand_cpg_gene_mappings(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands rows where 'unique_gene_name' contains multiple genes separated by commas.
    Also handles NaN values by labeling them as 'unmapped'.
    """
    expanded_df = df.copy()
    expanded_df['unique_gene_name'] = expanded_df['unique_gene_name'].fillna('unmapped')
    expanded_df['gene'] = expanded_df['unique_gene_name'].str.split(',')
    expanded_df = expanded_df.explode('gene')
    expanded_df['gene'] = expanded_df['gene'].str.strip()
    return expanded_df

def load_pi_cpg_mappings(file_path: str = "ewas_res_groupsig_128.xlsx") -> pd.DataFrame:
    """
    Loads the PI's provided CpG-to-gene mappings from Excel.
    """
    cols_to_load = ['cpg', 'chr', 'unique_gene_name', 'Start_hg38', 'End_hg38']
    df = pd.read_excel(file_path, usecols=cols_to_load)
    return expand_cpg_gene_mappings(df)

def attach_ewas_atlas_traits(mapping_df: pd.DataFrame, atlas_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregates trait associations from EWAS Atlas and attaches them to the mapping DataFrame.
    Implements rank scoring, correlation mapping, and strict formatting.
    """
    # 1. Map Correlations
    corr_map = {'pos': 'hyper', 'neg': 'hypo'}
    atlas_df = atlas_df.copy()
    atlas_df['correlation_mapped'] = atlas_df['correlation'].map(corr_map).fillna('NR')
    
    # 2. Calculate Rank Score (Numeric only if rank exists)
    def calc_rank_score(row):
        try:
            if pd.notna(row['rank']) and pd.notna(row['total_associations']) and row['total_associations'] != 0:
                score = float(row['rank']) / float(row['total_associations'])
                return score
        except:
            pass
        return None
    
    atlas_df['rank_score_num'] = atlas_df.apply(calc_rank_score, axis=1)
    
    # 3. Aggregation and Formatting
    def aggregate_traits(group):
        rows = []
        for _, row in group.iterrows():
            trait = str(row['trait']).strip()
            if not trait or trait.lower() == 'nan':
                continue
                
            # Get numeric score or default to 0.0 for formatting
            raw_score = row['rank_score_num']
            display_score = f"{raw_score:.3f}" if pd.notna(raw_score) else "0.000"
            
            corr = row['correlation_mapped']
            
            if pd.notna(row['pmid']):
                try:
                    pmid = str(int(float(row['pmid'])))
                except:
                    pmid = str(row['pmid'])
            else:
                pmid = "NA"
            
            rows.append({
                'trait': trait,
                'score_sort': raw_score if pd.notna(raw_score) else -1.0, # Put 0.000 hits at bottom
                'formatted': f"{trait}, {display_score}, {corr}, {pmid}"
            })
            
        # Sort by score descending, then trait alphabetical
        rows.sort(key=lambda x: (x['score_sort'], x['trait']), reverse=[True, False])
        
        trait_str = ";\n".join([r['formatted'] for r in rows]) if rows else None
        return pd.Series({'ewas_atlas_traits': trait_str})

    atlas_grouped = atlas_df.groupby('cpg').apply(aggregate_traits, include_groups=False).reset_index()
    result_df = pd.merge(mapping_df, atlas_grouped, on='cpg', how='left')
    return result_df

def analyze_unmapped_genes(atlas_df: pd.DataFrame, living_df: pd.DataFrame, original_mapping_df: pd.DataFrame):
    """
    Identifies genes in EWAS Atlas that are NOT in the living file symbols.
    Cross-checks with synonyms and separates established vs unestablished (decimals).
    """
    existing_symbols = set(living_df['symbol'].dropna().astype(str).unique())
    
    synonyms_map = {}
    for _, row in living_df.iterrows():
        symbol = row['symbol']
        synonyms_raw = row.get('synonyms')
        if pd.notna(synonyms_raw):
            try:
                if isinstance(synonyms_raw, str):
                    if synonyms_raw.startswith('['):
                        synonyms = json.loads(synonyms_raw)
                    else:
                        synonyms = [s.strip() for s in synonyms_raw.split(';') if s.strip()]
                else:
                    synonyms = synonyms_raw
                
                for syn in synonyms:
                    synonyms_map[str(syn)] = symbol
            except:
                pass

    original_data_map = original_mapping_df.set_index('cpg')['unique_gene_name'].to_dict()

    cpg_genes = atlas_df[['cpg', 'genes']].drop_duplicates()
    
    unmapped_records = []
    cpg_unmapped_genes = {}
    cpg_unmapped_regions = {}
    appendix_c_source = []
    
    for _, row in cpg_genes.iterrows():
        cpg = row['cpg']
        genes_str = row['genes']
        if pd.isna(genes_str):
            continue
            
        full_atlas_genes = [g.strip() for g in str(genes_str).split(';') if g.strip()]
        established = []
        unestablished = []
        
        for gene in full_atlas_genes:
            is_decimal = '.' in gene
            if is_decimal:
                appendix_c_source.append({'cpg': cpg, 'genes': gene})
            
            if gene in existing_symbols:
                continue
                
            parent_symbol = synonyms_map.get(gene)
            
            # Diagnostic Record
            orig_val = original_data_map.get(cpg, "")
            if isinstance(orig_val, str) and orig_val.lower() == 'unmapped':
                orig_val = ""
            elif pd.isna(orig_val):
                orig_val = ""

            unmapped_records.append({
                'cpg': cpg,
                'ewas_genes': genes_str,
                'uncaptured_gene': gene,
                'synonym_of': parent_symbol if parent_symbol else "",
                'nearest_from_original_data': orig_val,
                'is_decimal': is_decimal
            })
            
            if not parent_symbol:
                if is_decimal:
                    unestablished.append(gene)
                else:
                    established.append(gene)
        
        if established:
            cpg_unmapped_genes[cpg] = ";\n".join(sorted(established))
        if unestablished:
            cpg_unmapped_regions[cpg] = ";\n".join(sorted(unestablished))

    unaccounted_df = pd.DataFrame(unmapped_records)
    if not unaccounted_df.empty:
        unaccounted_df = unaccounted_df[unaccounted_df['is_decimal'] == False].drop(columns=['is_decimal'])
    else:
        unaccounted_df = pd.DataFrame(columns=['cpg', 'ewas_genes', 'uncaptured_gene', 'synonym_of', 'nearest_from_original_data'])
        
    appendix_c_df = pd.DataFrame(appendix_c_source).drop_duplicates()
    if appendix_c_df.empty:
        appendix_c_df = pd.DataFrame(columns=['cpg', 'genes'])

    return cpg_unmapped_genes, cpg_unmapped_regions, unaccounted_df, appendix_c_df
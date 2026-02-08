import pandas as pd
import numpy as np
import re
import json

def expand_cpg_gene_mappings(df: pd.DataFrame) -> pd.DataFrame:
    """
    Expands rows where 'unique_gene_name' contains multiple genes separated by commas.
    Also handles NaN values by labeling them as 'unmapped'.
    """
    # Create a copy to avoid modifying the original
    expanded_df = df.copy()
    
    # Fill NaN values in unique_gene_name
    expanded_df['unique_gene_name'] = expanded_df['unique_gene_name'].fillna('unmapped')
    
    # Convert unique_gene_name to list by splitting on comma
    expanded_df['gene'] = expanded_df['unique_gene_name'].str.split(',')
    
    # Explode the 'gene' list into separate rows
    expanded_df = expanded_df.explode('gene')
    
    # Clean up whitespace
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
    # pos → hyper; neg → hypo; NA → NR
    corr_map = {'pos': 'hyper', 'neg': 'hypo'}
    atlas_df = atlas_df.copy()
    atlas_df['correlation_mapped'] = atlas_df['correlation'].map(corr_map).fillna('NR')
    
    # 2. Calculate Rank Score
    # rank_score = rank / total_associations (only if rank exists)
    def calc_rank_score(row):
        try:
            if pd.notna(row['rank']) and pd.notna(row['total_associations']) and row['total_associations'] != 0:
                score = float(row['rank']) / float(row['total_associations'])
                return f"{score:.3f}"
        except:
            pass
        return "0.000"
    
    atlas_df['rank_score'] = atlas_df.apply(calc_rank_score, axis=1)
    
    # 3. Aggregation and Formatting
    def aggregate_traits(group):
        rows = []
        for _, row in group.iterrows():
            trait = str(row['trait']).strip()
            score = row['rank_score']
            corr = row['correlation_mapped']
            pmid = str(int(row['pmid'])) if pd.notna(row['pmid']) else "NA"
            
            rows.append({
                'trait': trait,
                'score': score,
                'formatted': f"{trait}, {score}, {corr}, {pmid}"
            })
            
        # Sort by rank_score descending, then trait alphabetical
        rows.sort(key=lambda x: (x['score'], x['trait']), reverse=[True, False])
        
        trait_str = ";\n".join([r['formatted'] for r in rows]) if rows else None
        
        return pd.Series({
            'ewas_atlas_traits': trait_str
        })

    atlas_grouped = atlas_df.groupby('cpg').apply(aggregate_traits, include_groups=False).reset_index()
    
    # Merge with mapping_df
    result_df = pd.merge(mapping_df, atlas_grouped, on='cpg', how='left')
    
    # No more n/a strings - leave as actual nulls/NaN
    
    return result_df

def analyze_unmapped_genes(atlas_df: pd.DataFrame, living_df: pd.DataFrame):
    """
    Identifies genes in EWAS Atlas that are NOT in the living file symbols.
    Cross-checks with synonyms and separates established vs unestablished (decimals).
    """
    # 1. Prepare set of symbols and synonyms map
    existing_symbols = set(living_df['symbol'].dropna().astype(str).unique())
    
    synonyms_map = {} # synonym -> original_symbol
    for _, row in living_df.iterrows():
        symbol = row['symbol']
        synonyms_raw = row.get('synonyms')
        if pd.notna(synonyms_raw):
            try:
                # Handle both JSON lists and string lists
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

    # 2. Extract genes per CpG from Atlas
    # ewas_atlas has 'cpg' and 'genes' (semicolon separated)
    cpg_genes = atlas_df[['cpg', 'genes']].drop_duplicates()
    
    unmapped_records = [] # For diagnostic table
    cpg_unmapped_strings = {} # For main XL column
    appendix_c_source = [] # For decimal genes
    
    for _, row in cpg_genes.iterrows():
        cpg = row['cpg']
        genes_str = row['genes']
        if pd.isna(genes_str):
            continue
            
        full_atlas_genes = [g.strip() for g in str(genes_str).split(';') if g.strip()]
        
        unmapped_established = []
        unmapped_unestablished = []
        
        for gene in full_atlas_genes:
            # Check if decimal (Unestablished)
            is_decimal = '.' in gene
            
            if is_decimal:
                appendix_c_source.append({'cpg': cpg, 'genes': gene})
                # We still process it for unmapped status? 
                # Spec: "exclude genes with decimals in them [from unaccounted_genes.csv]... create a separate table for those in appendix C"
                # Spec: "subset the genes that are not present anywhere in annotated_genes symbol... separate the list of genes by with-decimal"
            
            if gene in existing_symbols:
                continue
                
            # Check synonyms
            parent_symbol = synonyms_map.get(gene)
            
            # Diagnostic Record (Filtered decimals out later)
            unmapped_records.append({
                'cpg': cpg,
                'ewas_genes': genes_str,
                'uncaptured_gene': gene,
                'synonym_of': parent_symbol if parent_symbol else "None",
                'is_decimal': is_decimal
            })
            
            if not parent_symbol:
                # Truly unmapped (not symbol, not synonym)
                if is_decimal:
                    unmapped_unestablished.append(gene)
                else:
                    unmapped_established.append(gene)
        
        # Format the main XL column string
        if unmapped_established or unmapped_unestablished:
            parts = []
            if unmapped_established:
                parts.append(f"Established: {', '.join(sorted(unmapped_established))}")
            if unmapped_unestablished:
                parts.append(f"Unestablished: {', '.join(sorted(unmapped_unestablished))}")
            cpg_unmapped_strings[cpg] = "; ".join(parts)

    # 3. Finalize dataframes
    unaccounted_df = pd.DataFrame(unmapped_records)
    if not unaccounted_df.empty:
        # Filter OUT decimals from diagnostic table
        unaccounted_df = unaccounted_df[unaccounted_df['is_decimal'] == False].drop(columns=['is_decimal'])
    else:
        unaccounted_df = pd.DataFrame(columns=['cpg', 'ewas_genes', 'uncaptured_gene', 'synonym_of'])
        
    appendix_c_df = pd.DataFrame(appendix_c_source).drop_duplicates()
    if appendix_c_df.empty:
        appendix_c_df = pd.DataFrame(columns=['cpg', 'genes'])

    return cpg_unmapped_strings, unaccounted_df, appendix_c_df
import pandas as pd
import numpy as np

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
    """
    # Group atlas by CpG and aggregate traits and PMIDs
    def aggregate_traits(group):
        traits = sorted(list(set(group['trait'].dropna().astype(str))))
        pmids = sorted(list(set(group['pmid'].dropna().astype(int).astype(str))))
        
        trait_str = "; ".join(traits) if traits else "n/a"
        pmid_str = "; ".join(pmids) if pmids else "n/a"
        
        return pd.Series({
            'ewas_atlas_traits': trait_str,
            'ewas_atlas_pmids': pmid_str
        })

    atlas_grouped = atlas_df.groupby('cpg').apply(aggregate_traits, include_groups=False).reset_index()
    
    # Merge with mapping_df
    result_df = pd.merge(mapping_df, atlas_grouped, on='cpg', how='left')
    
    # Fill NaN for CpGs not in Atlas
    result_df['ewas_atlas_traits'] = result_df['ewas_atlas_traits'].fillna('n/a')
    result_df['ewas_atlas_pmids'] = result_df['ewas_atlas_pmids'].fillna('n/a')
    
    return result_df


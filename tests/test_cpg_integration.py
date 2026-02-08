import pytest
import pandas as pd
import sys
import os

# Add the project root to sys.path to import our modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.cpg_integration import expand_cpg_gene_mappings, load_pi_cpg_mappings, attach_ewas_atlas_traits

def test_cpg_mapping_expansion():
    # Sample data representing ewas_res_groupsig_128.xlsx
    data = {
        'cpg': ['cg07252486', 'cg22084806', 'cg06928226'],
        'chr': ['chr11', 'chr2', 'chr11'],
        'unique_gene_name': ['C11orf39, NTM', None, 'NTM'],
        'Start_hg38': [131663304, 113503692, 131644472],
        'End_hg38': [131663306, 113503694, 131644474]
    }
    df = pd.DataFrame(data)
    
    result = expand_cpg_gene_mappings(df)
    
    # Assertions
    assert len(result) == 4
    assert 'C11orf39' in result['gene'].values
    assert 'NTM' in result['gene'].values
    assert result[result['cpg'] == 'cg07252486'].shape[0] == 2
    assert result[result['cpg'] == 'cg22084806']['gene'].iloc[0] == 'unmapped'

def test_load_pi_cpg_mappings(tmp_path):
    # Create a temporary Excel file
    data = {
        'cpg': ['cg12345'],
        'chr': ['chr1'],
        'unique_gene_name': ['GENE1'],
        'Start_hg38': [100],
        'End_hg38': [200]
    }
    df = pd.DataFrame(data)
    file_path = tmp_path / "test_mappings.xlsx"
    df.to_excel(file_path, index=False)
    
    result = load_pi_cpg_mappings(str(file_path))
    
    assert len(result) == 1
    assert 'gene' in result.columns
    assert result['gene'].iloc[0] == 'GENE1'

def test_attach_ewas_atlas_traits():
    # Sample mapping data
    mapping_data = {
        'cpg': ['cg07252486', 'cg22084806', 'cg99999999'], # cg9999 is not in atlas
        'gene': ['NTM', 'unmapped', 'SOMEGENE']
    }
    mapping_df = pd.DataFrame(mapping_data)
    
    # Sample Atlas data
    atlas_data = {
        'cpg': ['cg07252486', 'cg07252486', 'cg22084806'],
        'trait': ['Body mass index', 'Obesity', 'Age'],
        'pmid': [12345, 67890, 11111]
    }
    atlas_df = pd.DataFrame(atlas_data)
    
    result = attach_ewas_atlas_traits(mapping_df, atlas_df)
    
    # Assertions
    assert 'ewas_atlas_traits' in result.columns
    # Traits for cg07252486 should be aggregated
    traits_cg1 = result[result['cpg'] == 'cg07252486']['ewas_atlas_traits'].iloc[0]
    assert 'Body mass index' in traits_cg1
    assert 'Obesity' in traits_cg1
    
    # Check cg22084806
    traits_cg2 = result[result['cpg'] == 'cg22084806']['ewas_atlas_traits'].iloc[0]
    assert 'Age' in traits_cg2

    # Check cg99999999 (not in atlas)
    traits_cg3 = result[result['cpg'] == 'cg99999999']['ewas_atlas_traits'].iloc[0]
    assert traits_cg3 == 'n/a'
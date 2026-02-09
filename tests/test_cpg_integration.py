import pytest
import pandas as pd
import numpy as np
import sys
import os

# Add the project root to sys.path to import our modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.cpg_integration import (
    expand_cpg_gene_mappings, 
    load_pi_cpg_mappings, 
    attach_ewas_atlas_traits,
    analyze_unmapped_genes
)

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

def test_attach_ewas_atlas_traits_refined():
    # Sample mapping data
    mapping_data = {
        'cpg': ['cg07252486', 'cg22084806', 'cg99999999'],
        'gene': ['NTM', 'unmapped', 'SOMEGENE']
    }
    mapping_df = pd.DataFrame(mapping_data)
    
    # Sample Atlas data
    atlas_data = {
        'cpg': ['cg07252486', 'cg07252486', 'cg22084806'],
        'trait': ['BMI', 'Obesity', 'Age'],
        'pmid': [123, 456, 789],
        'rank': [10, np.nan, 5],
        'total_associations': [100, 100, 50],
        'correlation': ['pos', 'neg', np.nan]
    }
    atlas_df = pd.DataFrame(atlas_data)
    
    result = attach_ewas_atlas_traits(mapping_df, atlas_df)
    
    # Assertions for formatting: trait, rank_score, correlation, pmid
    # cg07252486: 
    # BMI: score=0.100, corr=hyper, pmid=123
    # Obesity: score=0.000, corr=hypo, pmid=456
    traits_cg1 = result[result['cpg'] == 'cg07252486']['ewas_atlas_traits'].iloc[0]
    assert 'BMI, 0.100, hyper, 123' in traits_cg1
    assert 'Obesity, 0.000, hypo, 456' in traits_cg1
    
    # Check sorting: obesity (0.000) should be AFTER BMI (0.100) if descending score
    # Wait, spec says "ordered by rank_score, then subsequently alphabetical"
    # Actually, Obesity score 0.000, BMI score 0.100.
    # BMI should come first.
    lines = traits_cg1.split(';\n')
    assert 'BMI' in lines[0]
    
    # Check null handling: cg99999999 should be NaN/None, not "n/a"
    val = result[result['cpg'] == 'cg99999999']['ewas_atlas_traits'].iloc[0]
    assert pd.isna(val)

def test_analyze_unmapped_genes():
    # Living file data
    living_data = {
        'symbol': ['GENE1', 'GENE2'],
        'synonyms': ['["ALIAS1", "ALIAS2"]', None]
    }
    living_df = pd.DataFrame(living_data)
    
    # Original mapping data
    mapping_data = {
        'cpg': ['cg1', 'cg2'],
        'unique_gene_name': ['GENE1', 'ALIAS1']
    }
    mapping_df = pd.DataFrame(mapping_data)
    
    # Atlas data
    atlas_data = {
        'cpg': ['cg1', 'cg1', 'cg2', 'cg2', 'cg2'],
        'genes': ['GENE1;NEWGENE', 'GENE1;NEWGENE', 'ALIAS1;DECIMAL.1;ANOTHER', 'ALIAS1;DECIMAL.1;ANOTHER', 'ALIAS1;DECIMAL.1;ANOTHER']
    }
    atlas_df = pd.DataFrame(atlas_data)
    
    unmapped_genes, unmapped_regions, unaccounted_df, appendix_c_df = analyze_unmapped_genes(atlas_df, living_df, mapping_df)

    
    # unmapped_genes for cg1: NEWGENE
    assert unmapped_genes['cg1'] == 'NEWGENE'
    
    # unmapped_genes for cg2: ANOTHER
    # ALIAS1 should be found in synonyms of GENE1
    assert unmapped_genes['cg2'] == 'ANOTHER'
    
    # unmapped_regions for cg2: DECIMAL.1
    assert unmapped_regions['cg2'] == 'DECIMAL.1'
    
    # unaccounted_df: should have NEWGENE, ANOTHER, and ALIAS1
    assert 'NEWGENE' in unaccounted_df['uncaptured_gene'].values
    assert 'ALIAS1' in unaccounted_df['uncaptured_gene'].values
    assert 'DECIMAL.1' not in unaccounted_df['uncaptured_gene'].values
    
    # synonym_of for ALIAS1 should be GENE1
    row_alias = unaccounted_df[unaccounted_df['uncaptured_gene'] == 'ALIAS1']
    assert row_alias['synonym_of'].iloc[0] == 'GENE1'
    
    # appendix_c_df: should have DECIMAL.1
    assert 'DECIMAL.1' in appendix_c_df['genes'].iloc[0]


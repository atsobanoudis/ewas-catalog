import pandas as pd
import numpy as np
import pytest
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.disgenet_utils import annotate_with_disgenet

def test_annotate_with_disgenet_no_filter():
    # Sample gene data
    genes_df = pd.DataFrame({'symbol': ['GENE1', 'GENE2']})
    
    # Sample DisGeNET GEA data
    disgenet_data = {
        'gene_symbol': ['GENE1', 'GENE1', 'GENE2'],
        'disease_name': ['Disease A', 'Mental Disorder B', 'Disease C'],
        'diseaseClasses_UMLS_ST': ['Class 1', 'Mental or Behavioral Dysfunction (T048)', 'Class 2'],
        'score': [0.5, 0.7, 0.3],
        'polarity': ['Positive', 'Positive', 'Positive'],
        'pmYear': [2020, 2021, 2019],
        'reference_type': ['PMID', 'PMID', 'PMID'],
        'reference': [123, 456, 789],
        'source': ['S1', 'S2', 'S3'],
        'associationType': ['T1', 'T2', 'T3']
    }
    disgenet_df = pd.DataFrame(disgenet_data)
    
    # Run WITHOUT psych filter
    result = annotate_with_disgenet(genes_df, disgenet_df, psych_only=False)
    
    # GENE1 should have BOTH diseases
    diseases_gene1 = result[result['symbol'] == 'GENE1']['disgenet_diseases'].iloc[0]
    assert 'Mental Disorder B, 0.7' in diseases_gene1
    assert 'Disease A, 0.5' in diseases_gene1
    
    # GENE2 should have its disease
    diseases_gene2 = result[result['symbol'] == 'GENE2']['disgenet_diseases'].iloc[0]
    assert 'Disease C, 0.3' in diseases_gene2

def test_annotate_with_disgenet_psych_filter():
    # Sample gene data
    genes_df = pd.DataFrame({'symbol': ['GENE1', 'GENE2']})
    
    # Sample DisGeNET GEA data
    disgenet_data = {
        'gene_symbol': ['GENE1', 'GENE1', 'GENE2'],
        'disease_name': ['Disease A', 'Mental Disorder B', 'Disease C'],
        'diseaseClasses_UMLS_ST': ['Class 1', 'Mental or Behavioral Dysfunction (T048)', 'Class 2'],
        'score': [0.5, 0.7, 0.3],
        'polarity': ['Positive', 'Positive', 'Positive'],
        'pmYear': [2020, 2021, 2019],
        'reference_type': ['PMID', 'PMID', 'PMID'],
        'reference': [123, 456, 789],
        'source': ['S1', 'S2', 'S3'],
        'associationType': ['T1', 'T2', 'T3']
    }
    disgenet_df = pd.DataFrame(disgenet_data)
    
    # Run WITH psych filter
    result = annotate_with_disgenet(genes_df, disgenet_df, psych_only=True)
    
    # GENE1 should only have Mental Disorder B
    diseases_gene1 = result[result['symbol'] == 'GENE1']['disgenet_psych_diseases'].iloc[0]
    assert 'Mental Disorder B, 0.7' in diseases_gene1
    assert 'Disease A' not in diseases_gene1
    
    # GENE2 should be NaN or 'n/a'
    diseases_gene2 = result[result['symbol'] == 'GENE2']['disgenet_psych_diseases'].iloc[0]
    assert pd.isna(diseases_gene2) or diseases_gene2 == 'n/a'

def test_annotate_with_disgenet_edge_cases():
    # 1. Missing diseaseClasses_UMLS_ST column
    genes_df = pd.DataFrame({'symbol': ['GENE1', 'GENE_NOT_FOUND']})
    disgenet_df = pd.DataFrame({
        'gene_symbol': ['GENE1'],
        'disease_name': ['D1'],
        'score': [0.5],
        'polarity': ['Positive'],
        'pmYear': [2020],
        'reference_type': ['PMID'],
        'reference': [123],
        'source': ['S1'],
        'associationType': ['T1']
    })
    result = annotate_with_disgenet(genes_df, disgenet_df, psych_only=True)
    assert 'disgenet_psych_diseases' in result.columns
    # GENE_NOT_FOUND should hit the empty branch (lines 44-46)
    assert pd.isna(result[result['symbol'] == 'GENE_NOT_FOUND']['disgenet_psych_diseases'].iloc[0])

    # 2. Multi-score error
    disgenet_df_multi = pd.DataFrame({
        'gene_symbol': ['GENE1', 'GENE1'],
        'disease_name': ['D1', 'D1'],
        'score': [0.5, 0.6], # Different scores for same disease
        'polarity': ['Positive', 'Positive'],
        'pmYear': [2020, 2021],
        'reference_type': ['PMID', 'PMID'],
        'reference': [123, 456],
        'source': ['S1', 'S1'],
        'associationType': ['T1', 'T1']
    })
    result_multi = annotate_with_disgenet(genes_df, disgenet_df_multi, psych_only=False)
    assert 'D1, error' in result_multi['disgenet_diseases'].iloc[0]

    # 3. Non-PMID reference and NA reference
    disgenet_df_refs = pd.DataFrame({
        'gene_symbol': ['GENE1', 'GENE1'],
        'disease_name': ['D1', 'D2'],
        'score': [0.5, 0.5],
        'polarity': [np.nan, 'Positive'],
        'pmYear': [np.nan, 2020],
        'reference_type': ['Other', 'PMID'],
        'reference': [np.nan, 'abc'], # Testing non-int PMID and NA other
        'source': ['S1', 'S2'],
        'associationType': ['T1', 'T2']
    })
    result_refs = annotate_with_disgenet(genes_df, disgenet_df_refs, psych_only=False)
    evidence = result_refs['disgenet_evidence'].iloc[0]
    assert 'NA (S1, T1)' in evidence
    assert 'abc' in evidence
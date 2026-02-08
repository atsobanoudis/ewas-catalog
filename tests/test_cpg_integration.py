import pytest
import pandas as pd
import sys
import os

# Add the project root to sys.path to import our modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.cpg_integration import expand_cpg_gene_mappings, load_pi_cpg_mappings



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

        'End_hg38': [200],

        'extra_col': ['ignore_me']

    }

    df = pd.DataFrame(data)

    file_path = tmp_path / "test_mappings.xlsx"

    df.to_excel(file_path, index=False)

    

    result = load_pi_cpg_mappings(str(file_path))

    

    assert len(result) == 1

    assert 'gene' in result.columns

    assert result['gene'].iloc[0] == 'GENE1'

    assert 'extra_col' not in result.columns



import pandas as pd
import pytest
import sys
import os
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.viz_engine import plot_disease_heatmap

def test_plot_disease_heatmap(tmp_path):
    # Sample data
    data = {
        'symbol': ['GENE1', 'GENE1', 'GENE2', 'GENE3'],
        'disgenet_diseases': [
            'Schizophrenia, 0.7;\nBipolar Disorder, 0.5',
            'Schizophrenia, 0.7;\nBipolar Disorder, 0.5',
            'Schizophrenia, 0.4',
            'Depression, 0.3'
        ]
    }
    df = pd.DataFrame(data)
    
    output_file = tmp_path / "test_heatmap.png"
    
    # This should run without error and create a file
    plot_disease_heatmap(df, output_path=str(output_file))
    
    assert Path(output_file).is_file()
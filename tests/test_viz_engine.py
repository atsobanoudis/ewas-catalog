import pandas as pd
import pytest
import sys
import os
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.viz_engine import plot_disease_heatmap, plot_broad_spectrum_network

def test_plot_disease_heatmap(tmp_path):
    # ... (existing test)
    pass

def test_plot_broad_spectrum_network(tmp_path):
    # Sample data
    data = {
        'symbol': ['GENE1', 'GENE2'],
        'disgenet_diseases': [
            'Disease A, 0.7;\nDisease B, 0.5',
            'Disease A, 0.4'
        ]
    }
    df = pd.DataFrame(data)
    
    output_file = tmp_path / "test_network.png"
    
    # This should run without error and create a file
    plot_broad_spectrum_network(df, output_path=str(output_file))
    
    assert Path(output_file).is_file()
from original_annotation.modules.viz_engine import plot_broad_spectrum_network
import pandas as pd
import pytest
from pathlib import Path
import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def test_network_refinements(tmp_path):
    # Data with varying scores and long names
    data = {
        'symbol': ['GENE1', 'GENE2'],
        'disgenet_diseases': [
            'Very Long Disease Name That Should Be Wrapped, 0.05;\nHigh Score Disease, 0.8',
            'Low Score Disease, 0.01'
        ]
    }
    df = pd.DataFrame(data)
    
    output_file = tmp_path / "test_refined_network.png"
    
    # Test with filtering threshold 0.1
    # 'Very Long...' (0.05) and 'Low Score...' (0.01) should be filtered out
    # 'High Score Disease' (0.8) should remain
    # Text wrapping should be applied to 'High Score Disease' (though it's short, logic runs)
    
    plot_broad_spectrum_network(df, output_path=str(output_file), score_threshold=0.1)
    
    assert Path(output_file).is_file()
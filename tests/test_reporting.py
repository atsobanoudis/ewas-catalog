import pandas as pd
import pytest
import sys
import os
from pathlib import Path

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from original_annotation.modules.reporting_engine import generate_summary_report

def test_generate_summary_report(tmp_path):
    # Sample augmented data
    data = {
        'symbol': ['GENE1', 'GENE2', 'GENE3'],
        'disgenet_psych_diseases': [
            'Schizophrenia, 0.7;\nBipolar Disorder, 0.5',
            'Depression, 0.2',
            'n/a'
        ],
        'ewas_atlas_traits': [
            'BMI;\nSmoking',
            'n/a',
            'Age'
        ]
    }
    df = pd.DataFrame(data)
    
    output_file = tmp_path / "summary_report.md"
    
    # Generate report
    generate_summary_report(df, output_path=str(output_file))
    
    assert Path(output_file).is_file()
    
    # Check content
    content = Path(output_file).read_text()
    assert "# Psychiatric Epigenetics Discovery Report" in content
    assert "Top Psychiatric Associations" in content
    assert "GENE1" in content
    assert "Schizophrenia" in content
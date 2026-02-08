# Tech Stack

## Programming Languages
- **Python**: Primary language for data aggregation, API interactions, and pipeline orchestration.
- **R**: Used for specific statistical analysis, DisGeNET data pulls, and environment-specific tasks.

## Data Processing & Analysis
- **pandas (Python)**: Core library for data manipulation, cleaning, and Excel/CSV handling.
- **requests (Python)**: Used for all REST API interactions (NCBI, HGNC, UniProt, etc.).
- **tidyverse (R)**: (Inferred) Likely used for data transformation and visualization in R.

## Environment & Dependency Management
- **renv (R)**: Manages R package versions and ensures reproducibility.
- **venv (Python)**: Standard Python virtual environment for managing project-specific packages.

## Data Storage & Formats
- **Excel (.xlsx)**: Primary output format for annotated gene lists (`annotated_genes.xlsx`).
- **CSV/TSV**: Used for intermediate datasets and database exports (e.g., GWAS catalog).
- **RDS (R)**: Used for storing structured R data objects (e.g., DisGeNET associations).
- **JSON**: Used for local data caching and configuration.

## Reporting & Visualization
- **Markdown**: For project documentation and structured progress reports.
- **Matplotlib/Seaborn**: For generating publication-ready static heatmaps (`viz_engine.py`).
- **NetworkX**: For generating gene-disease network graphs (`viz_engine.py`).
- **disgenet2r (R)**: For specialized DisGeNET visualizations.


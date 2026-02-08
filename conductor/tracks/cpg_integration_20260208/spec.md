# Specification: Integrate CpG-Epigenetic Data and Implement Core Visualization Suite

## Overview
This track aims to finalize the integration of the PI's raw CpG data (chromosomal locations and nearest genes) into the existing `annotated_genes.xlsx` pipeline. Additionally, it will introduce a core visualization module to generate publication-quality heatmaps and networks, satisfying the goal of providing "cerebral" interpretation and "wow-factor" results for the PI.

## Objectives
1.  **Data Integration:**
    -   Update the annotation pipeline to include specific CpG chromosomal locations.
    -   Integrate CpG-to-trait associations (all traits) from `ewas_atlas.csv`.
    -   Expand gene-disease associations to include non-mental health traits by running a non-filtered version of the DisGeNET pipeline.
    -   Ensure a clean mapping between CpG sites and their annotated gene metadata.
2.  **Visualization Module:**
    -   Develop a Python module for generating "Disease Overlap Heatmaps" using data from GWAS, DisGeNET, and PubMed.
    -   Implement basic aesthetic customization as per Product Guidelines.
3.  **Reporting:**
    -   Update `REPORT.md` (or generate a new structured report) that summarizes these findings with an interpretive narrative.

## Scope
-   In-scope: CpG site mapping update, EWAS Atlas trait integration, Full-spectrum DisGeNET association appending, Heatmap generation, Summary reporting.
-   Out-of-scope: Interactive dashboard (web-based), Full manuscript drafting.

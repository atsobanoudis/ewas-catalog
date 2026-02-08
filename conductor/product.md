# Initial Concept

The project originated as a request to annotate a list of CpG sites (identified via chromosomal locations) with their nearest genes. The goal was to transform raw epigenetic data into a structured `annotated_genes.xlsx` file. This has since evolved into a more ambitious "discovery suite" that integrates multi-omic data (gene metadata, functional summaries, and disease associations) to facilitate high-level interpretation and publication-ready reporting for psychiatric epigenetics research.

# Product Guide

## Vision
To elevate a standard gene annotation task into a comprehensive, high-impact **Psychiatric Epigenetics Discovery Suite**. This tool not only aggregates data but synthesizes it into "publication-ready" insights, directly supporting hypothesis generation and manuscript preparation. The ultimate goal is to exceed the PI's expectations by delivering professional-grade data visualizations and automated reports that secure co-authorship by demonstrating significant intellectual contribution.

## Target Audience
- **Primary:** Principal Investigator (PI).
- **Secondary:** The User (Medical Student Researcher) - focusing on efficiency, reproducibility, and "wow-factor" results.
- **Tertiary:** Future lab members or collaborators who need to replicate or extend the analysis.

## Core Objectives
1.  **Robust Annotation Pipeline:** Maintain and enhance the current `annotated_genes.xlsx` generation, ensuring it is reliable and incorporates the newly requested epigenetic data (CpG locations, nearby genes).
2.  **Publication-Ready Visualization:** Implement high-impact visualizations (e.g., Disease Overlap Heatmaps, Functional Enrichment Networks) that provide immediate "cerebral" value and interpretation of the data, requiring minimal effort from the PI to understand.
3.  **Automated Insight Reporting:** Generate polished reports (PDF/HTML) that synthesize findings, providing a narrative layer on top of the raw data to satisfy ICMJE authorship criteria for intellectual contribution.
4.  **Epigenetic Context:** Seamlessly integrate the provided CpG data, visualizing relationships between specific CpG sites and gene function/disease association.

## Key Features
-   **Multi-Source Data Aggregation:** Automated fetching and cleaning of data from NCBI, HGNC, UniProt, GWAS Catalog, Harmonizome, PubMed, and DisGeNET.
-   **Intellectual Synthesis Layer:** Algorithms to score or rank gene-disease associations based on evidence strength, moving beyond simple listing.
-   **"One-Click" Reporting:** A command that outputs a formatted summary of top findings, suitable for immediate insertion into a grant or paper.
-   **Interactive/Static Visuals:** Generation of heatmaps and network graphs that highlight clustering of genes around specific psychiatric traits.

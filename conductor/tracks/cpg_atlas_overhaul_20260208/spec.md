# Specification: Enhanced CpG Integration, EWAS Atlas Refinement, and Report Overhaul

## Overview
This track upgrades the data augmentation process to ensure total data integrity, refined EWAS Atlas integration, and comprehensive gap analysis. It also overhauls `REPORT.md` to document these new methodologies and present critical findings.

## Functional Requirements

### 1. Data Processing Refinement (`augment_living_file.py` / `cpg_integration.py`)
-   **Comprehensive Join:** Ensure *every* unique CpG from `ewas_res_groupsig_128.xlsx` appears in the output.
    -   CpGs with no gene mapping must have blank (null) cells for genetic columns.
-   **EWAS Atlas Enhancement:**
    -   **Rank Score:** Calculate `rank_score = rank / total_associations` (3 decimal places). If rank is missing, `rank_score = 0.000`.
    -   **Correlation Mapping:** `pos` → `hyper`, `neg` → `hypo`, `NA` → `NR`.
    -   **Formatting:** `trait, rank_score, correlation, pmid`.
    -   **Sorting:** Descending by `rank_score`, then alphabetical.
-   **Unmapped Genes Analysis:**
    -   **Logic:** For each CpG, identify genes in `ewas_atlas.csv` not present in `annotated_genes.xlsx` `symbol`.
    -   **Synonym Check:** Cross-reference these potential new genes against the `synonyms` column.
    -   **`ewas_unmapped_gene` Column (Main XL):**
        -   Contains genes that are **NOT** in `symbol` **AND NOT** in `synonyms`.
        -   Format: `Established: GeneA; Unestablished: GeneB.1` (Decimal separation).
    -   **Unaccounted Genes Table (Diagnostic CSV):**
        -   Rows for every gene in Atlas not in `symbol` (filtered: No Decimals).
        -   Columns: `cpg`, `ewas_genes` (full list), `uncaptured_gene` (the specific one), `synonym_of` (Symbol if found in synonyms, else None).
    -   **Appendix C Source (Diagnostic CSV):**
        -   Rows for every gene in Atlas not in `symbol` (filtered: Decimals Only).

### 2. Documentation Overhaul (`REPORT.md`)
-   **Detailed Workflow (EWAS Atlas):**
    -   Update to include `rank_score`, formatting, and the unmapped gene logic (including synonym check).
-   **Interesting Findings:**
    -   **FMN1 Case Study:** Insert "FMN1" discussion with specific genomic coordinates and trait details.
    -   **Unaccounted Genes Table:** Insert/Reference the generated table of uncaptured genes (and their synonym status).
-   **Considerations/Limitations:**
    -   Discuss methodology (rank scores, p-values).
    -   Discuss "Nearest Gene" vs "EWAS-Associated Gene".
    -   Use `[[INSERT_DETAIL: ...]]` for manual review items.
-   **Appendix C:**
    -   Insert table of decimal-named genes.

## Acceptance Criteria
-   `annotated_genes_augmented.xlsx` has `ewas_unmapped_gene` column (verified vs synonyms).
-   `unaccounted_genes.csv` is generated with `synonym_of` column (no decimals).
-   `appendix_c_source.csv` is generated (decimals only).
-   `REPORT.md` updated with placeholders.

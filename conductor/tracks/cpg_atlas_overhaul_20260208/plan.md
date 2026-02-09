# Implementation Plan: Enhanced CpG Integration & Report Overhaul

## Phase 1: Data Augmentation Logic [checkpoint: b688562]
- [x] Task: Update `original_annotation/modules/cpg_integration.py` to implement refined EWAS Atlas logic. [b688562]
    - [x] Implement `rank_score` calculation and correlation mapping.
    - [x] Implement robust synonym checking for unmapped genes.
    - [x] Separate `ewas_unmapped_gene` (Established) and `ewas_unmapped_regions` (Unestablished/Decimal).
    - [x] Remove prefixes from unmapped column strings and handle empty associations as nulls.
    - [x] Create diagnostic tables (Unaccounted Genes vs Appendix C Source).
- [x] Task: Update `augment_living_file.py` to use the comprehensive join and new Atlas features. [b688562]
    - [x] Switch to a comprehensive join (Outer/Right) to preserve all master CpGs.
    - [x] Integrate the new `ewas_unmapped_gene` and `ewas_unmapped_regions` columns.
    - [x] Fix duplicate `disgenet_diseases` column issue.
    - [x] Output `unaccounted_genes.csv` (Filtered: No Decimals) and `appendix_c_source.csv` (Filtered: Decimals Only).
- [x] Task: Conductor - User Manual Verification 'Phase 1: Data Augmentation Logic' (Protocol in workflow.md)

## Phase 2: Documentation Overhaul [checkpoint: 3721ef4]
- [x] Task: Update `REPORT.md` - Detailed Workflow. [3721ef4]
    - [x] Add `rank_score`, correlation mapping, and unmapped gene logic to the EWAS Atlas section.
- [x] Task: Update `REPORT.md` - Interesting Findings. [3721ef4]
    - [x] Insert the **FMN1 Case Study** with placeholders for specific stats.
    - [x] Insert the **Unaccounted Genes** findings (referencing the CSV).
- [x] Task: Update `REPORT.md` - Considerations/Limitations. [3721ef4]
    - [x] Discuss methodology (rank limitations, p-value exclusion).
    - [x] Discuss scientific implications (Nearest vs. Associated).
- [x] Task: Update `REPORT.md` - Appendix C. [3721ef4]
    - [x] Create Appendix C table for decimal-named genes (Unestablished).
- [x] Task: Conductor - User Manual Verification 'Phase 2: Documentation Overhaul' (Protocol in workflow.md)

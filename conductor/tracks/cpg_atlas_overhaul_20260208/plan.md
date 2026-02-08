# Implementation Plan: Enhanced CpG Integration & Report Overhaul

## Phase 1: Data Augmentation Logic [checkpoint: 7f89e3c]
- [x] Task: Update `original_annotation/modules/cpg_integration.py` to implement refined EWAS Atlas logic. [7f89e3c]
    - [x] Implement `rank_score` calculation and correlation mapping.
    - [x] Implement robust synonym checking for unmapped genes.
    - [x] Create function to generate `ewas_unmapped_gene` string (Established/Unestablished split).
    - [x] Create function to generate diagnostic tables (Unaccounted Genes vs Appendix C Source).
- [x] Task: Update `augment_living_file.py` to use the comprehensive join and new Atlas features. [7f89e3c]
    - [x] Switch to a comprehensive join (Outer/Right) to preserve all master CpGs.
    - [x] Integrate the new `ewas_unmapped_gene` column generation.
    - [x] Output `unaccounted_genes.csv` (Filtered: No Decimals).
    - [x] Output `appendix_c_source.csv` (Filtered: Decimals Only).
- [x] Task: Conductor - User Manual Verification 'Phase 1: Data Augmentation Logic' (Protocol in workflow.md)

## Phase 2: Documentation Overhaul
- [ ] Task: Update `REPORT.md` - Detailed Workflow.
    - [ ] Add `rank_score`, correlation mapping, and unmapped gene logic to the EWAS Atlas section.
- [ ] Task: Update `REPORT.md` - Interesting Findings.
    - [ ] Insert the **FMN1 Case Study** with placeholders for specific stats.
    - [ ] Insert the **Unaccounted Genes** findings (referencing the CSV).
- [ ] Task: Update `REPORT.md` - Considerations/Limitations.
    - [ ] Discuss methodology (rank limitations, p-value exclusion).
    - [ ] Discuss scientific implications (Nearest vs. Associated).
- [ ] Task: Update `REPORT.md` - Appendix C.
    - [ ] Create Appendix C table for decimal-named genes (Unestablished).
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Documentation Overhaul' (Protocol in workflow.md)

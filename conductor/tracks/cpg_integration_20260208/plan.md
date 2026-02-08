# Implementation Plan: Integrate CpG-Epigenetic Data and Implement Core Visualization Suite

## Phase 1: Data Integration & Broad Association Mapping
- [x] Task: Integrate PI-provided CpG-to-gene mappings from `ewas_res_groupsig_128.xlsx`. [f3d50b4]
    - [x] Ensure all CpG chromosomal coordinates and their pre-mapped genes are correctly loaded.
- [x] Task: Integrate CpG-to-trait associations from `ewas_atlas.csv`. [9cd7345]
    - [x] Match by CpG ID; for gene-level traits, defer to the PI's provided mappings for consistency.
- [x] Task: Implement broad-spectrum gene-disease associations. [4fa3864]
    - [x] Update `disgenet_append_gea.py` (or create a variant) to include associations without filtering for "Mental or Behavioral Dysfunction".
- [ ] Task: Modify `original-annotation/main.py` or equivalent to include CpG locations and broad trait associations in output.
- [ ] Task: Conductor - User Manual Verification 'Phase 1: Data Integration & Broad Association Mapping' (Protocol in workflow.md)

## Phase 2: Core Visualization Suite
- [ ] Task: Design and implement the Visualization Module (`viz_engine.py`).
    - [ ] Implement `plot_disease_heatmap(gene_list, association_data)` using Seaborn/Matplotlib.
    - [ ] Add aesthetic configuration support (palettes, fonts) to match Academic & Clean guidelines.
- [ ] Task: Generate initial heatmaps for the current `annotated_genes.xlsx` dataset.
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Core Visualization Suite' (Protocol in workflow.md)

## Phase 3: Automated Reporting & Synthesis
- [ ] Task: Create a report generation script to synthesize findings.
    - [ ] Implement logic to identify "top-hit" gene-disease clusters for interpretive summaries.
    - [ ] Output a formatted Markdown or HTML summary.
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Automated Reporting & Synthesis' (Protocol in workflow.md)

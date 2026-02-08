# Implementation Plan: Integrate CpG-Epigenetic Data and Implement Core Visualization Suite

## Phase 1: Data Integration & Broad Association Mapping [checkpoint: db91993]
- [x] Task: Integrate PI-provided CpG-to-gene mappings from `ewas_res_groupsig_128.xlsx`. [f3d50b4]
    - [x] Ensure all CpG chromosomal coordinates and their pre-mapped genes are correctly loaded.
- [x] Task: Integrate CpG-to-trait associations from `ewas_atlas.csv`. [9cd7345]
    - [x] Match by CpG ID; for gene-level traits, defer to the PI's provided mappings for consistency.
- [x] Task: Implement broad-spectrum gene-disease associations. [4fa3864]
    - [x] Update `disgenet_append_gea.py` (or create a variant) to include associations without filtering for "Mental or Behavioral Dysfunction".
- [x] Task: Create `augment_living_file.py` to update the existing `annotated_genes.xlsx` with CpG info and broad associations. [ca4d92f]
- [x] Task: Conductor - User Manual Verification 'Phase 1: Data Integration & Broad Association Mapping' (Protocol in workflow.md)


## Phase 2: Core Visualization Suite
- [x] Task: Design and implement the Visualization Module (`viz_engine.py`). [566a328]
    - [x] Implement `plot_disease_heatmap(gene_list, association_data)` using Seaborn/Matplotlib.
    - [x] Add aesthetic configuration support (palettes, fonts) to match Academic & Clean guidelines.
- [x] Task: Generate initial heatmaps for the current `annotated_genes.xlsx` dataset. [c9f99f0]

## Phase 2.1: Visualization Refinements (Network & R Integration)
- [x] Task: Update `viz_engine.py` to revert psychiatric heatmap to 300 DPI. [f531551]
    - [x] Add `dpi` parameter to `plot_disease_heatmap`.
- [x] Task: Implement `plot_broad_spectrum_network` in `viz_engine.py` using `networkx` to replace the broad-spectrum heatmap. [474965c]
- [x] Task: Create/Update an R script (`disgenet/plot_heatmap.R`) to generate a static 300 DPI heatmap using `disgenet2r` and the existing RDS file. [a041772]
- [x] Task: Refine network graph (weights, text wrapping, filtering) and fix R script printing. [cf6e2d1]
- [x] Task: Integrate these new visualizations into `augment_living_file.py`. [b823429]
- [ ] Task: Conductor - User Manual Verification 'Phase 2.1: Visualization Refinements' (Protocol in workflow.md)
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Core Visualization Suite' (Protocol in workflow.md)

## Phase 3: Automated Reporting & Synthesis
- [ ] Task: Create a report generation script to synthesize findings.
    - [ ] Implement logic to identify "top-hit" gene-disease clusters for interpretive summaries.
    - [ ] Output a formatted Markdown or HTML summary.
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Automated Reporting & Synthesis' (Protocol in workflow.md)

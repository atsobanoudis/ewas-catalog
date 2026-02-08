import pandas as pd
from pathlib import Path
import sys
import os
import subprocess
import shutil

# Add root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from original_annotation.modules.cpg_integration import load_pi_cpg_mappings, attach_ewas_atlas_traits, analyze_unmapped_genes
from original_annotation.modules.disgenet_utils import annotate_with_disgenet
from original_annotation.modules.viz_engine import plot_disease_heatmap, plot_broad_spectrum_network
from original_annotation.modules.reporting_engine import generate_summary_report


def augment_living_file(living_path: str = "annotated_genes.xlsx", 
                        mapping_path: str = "ewas_res_groupsig_128.xlsx",
                        atlas_path: str = "data/ewas_atlas.csv",
                        gea_path: str = "data/disgenet_gea.csv",
                        output_path: str = "annotated_genes_augmented.xlsx"):
    """
    Loads the existing living file and appends CpG coordinates, Atlas traits, and broad DisGeNET data.
    Ensures all CpGs are preserved and performs gap analysis on unmapped genes.
    """
    print(f"[AUGMENT] Loading living file: {living_path}")
    df_living = pd.read_excel(living_path)
    
    print(f"[AUGMENT] Loading CpG mappings from PI file: {mapping_path}")
    df_mappings = load_pi_cpg_mappings(mapping_path)
    
    # 1. Join genetic annotations to the CpG mappings (Comprehensive Join)
    # We use a right join on mappings to ensure every CpG provided by PI is retained
    print("[AUGMENT] Performing comprehensive join (retaining all CpGs)...")
    # Clean mappings to just what we need for join
    df_mappings_clean = df_mappings[['cpg', 'chr', 'Start_hg38', 'End_hg38', 'gene']].rename(columns={
        'chr': 'cpg_chr',
        'Start_hg38': 'cpg_start',
        'End_hg38': 'cpg_end'
    })
    
    # Ensure 'input' is clean
    df_living['input'] = df_living['input'].astype(str).str.strip()
    
    # Right join: keep all rows from mappings_clean, even if no match in living
    df_augmented = pd.merge(df_living, df_mappings_clean, left_on='input', right_on='gene', how='right')
    
    # Drop the redundant 'gene' column
    if 'gene' in df_augmented.columns:
        df_augmented = df_augmented.drop(columns=['gene'])

    # 2. EWAS Atlas Integration & Gap Analysis
    if Path(atlas_path).is_file():
        print("[AUGMENT] Processing EWAS Atlas data and gap analysis...")
        df_atlas = pd.read_csv(atlas_path, low_memory=False)
        
        # A. Attach traits (with refined formatting/scoring)
        df_augmented = attach_ewas_atlas_traits(df_augmented, df_atlas)
        
        # B. Analyze unmapped genes
        unmapped_map, unaccounted_df, appendix_c_df = analyze_unmapped_genes(df_atlas, df_living)
        
        # C. Attach ewas_unmapped_gene column
        df_augmented['ewas_unmapped_gene'] = df_augmented['cpg'].map(unmapped_map)
        
        # D. Export diagnostic tables
        print("[AUGMENT] Exporting gap analysis diagnostic tables...")
        unaccounted_df.to_csv("unaccounted_genes.csv", index=False)
        appendix_c_df.to_csv("appendix_c_source.csv", index=False)
    else:
        print(f"[AUGMENT] Warning: {atlas_path} not found.")

    # 3. Rename existing 'disgenet_evidence' to 'disgenet_psych_evidence' if it exists
    if 'disgenet_evidence' in df_augmented.columns:
        print("[AUGMENT] Renaming 'disgenet_evidence' to 'disgenet_psych_evidence'...")
        df_augmented = df_augmented.rename(columns={'disgenet_evidence': 'disgenet_psych_evidence'})

    # 4. Attach Broad-Spectrum DisGeNET
    if Path(gea_path).is_file():
        print("[AUGMENT] Attaching Broad-Spectrum DisGeNET associations...")
        df_gea = pd.read_csv(gea_path)
        # Note: annotate_with_disgenet handles both psych and broad spectrum.
        # We target 'symbol' column. Rows without genes (NaN symbols) will be skipped by the loop in disgenet_utils.
        df_augmented = annotate_with_disgenet(df_augmented, df_gea, psych_only=False)
        
        # Clean up unwanted evidence column for broad spectrum
        if 'disgenet_evidence' in df_augmented.columns:
             df_augmented = df_augmented.drop(columns=['disgenet_evidence'])
    else:
        print(f"[AUGMENT] Warning: {gea_path} not found.")

    # Reorder columns: CpG info first, then genetic info, then Atlas traits, then unmapped genes
    cols = list(df_augmented.columns)
    cpg_cols = ['cpg', 'cpg_chr', 'cpg_start', 'cpg_end']
    for c in reversed(cpg_cols):
        if c in cols:
            cols.remove(c)
            cols.insert(0, c)
    
    # Ensure unmapped gene column is near the end
    if 'ewas_unmapped_gene' in cols:
        cols.remove('ewas_unmapped_gene')
        cols.append('ewas_unmapped_gene')
        
    df_augmented = df_augmented[cols]

    print(f"[AUGMENT] Saving augmented file to: {output_path}")
    df_augmented.to_excel(output_path, index=False)


    # 4. Generate Visualizations
    print("[AUGMENT] Generating visualizations...")
    viz_dir = Path("viz")
    viz_dir.mkdir(exist_ok=True)
    
    # Psychiatric Heatmap (Python - Seaborn) - Reverted to 300 DPI
    if 'disgenet_psych_diseases' in df_augmented.columns:
        plot_disease_heatmap(df_augmented, 
                             disease_col='disgenet_psych_diseases', 
                             output_path="viz/heatmap_psychiatric_seaborn.png",
                             title="Psychiatric Gene-Disease Associations (Seaborn)",
                             dpi=300)
    
    # Broad-Spectrum Network Graph (Python - NetworkX)
    if 'disgenet_diseases' in df_augmented.columns:
        plot_broad_spectrum_network(df_augmented, 
                                    disease_col='disgenet_diseases', 
                                    output_path="viz/network_broad_spectrum.png",
                                    title="Broad-Spectrum Gene-Disease Network",
                                    score_threshold=0.1, # Filter low scores
                                    wrap_width=20) # Wrap long labels

    # Psychiatric Heatmap (R - DisGeNET2R)
    # Check if Rscript is available
    if shutil.which("Rscript"):
        r_script_path = Path("disgenet/plot_heatmap.R")
        if r_script_path.is_file():
            print("[AUGMENT] Executing R script for DisGeNET2R heatmap...")
            try:
                subprocess.run(["Rscript", str(r_script_path)], check=True)
            except subprocess.CalledProcessError as e:
                print(f"[AUGMENT] Error running R script: {e}")
        else:
            print(f"[AUGMENT] R script not found at {r_script_path}")
    else:
        print("[AUGMENT] Rscript not found. Skipping DisGeNET2R heatmap generation.")

    # 5. Generate Summary Report
    print("[AUGMENT] Generating summary report...")
    generate_summary_report(df_augmented, output_path="summary_report.md")

    print("[AUGMENT] Done.")

if __name__ == "__main__":
    augment_living_file()
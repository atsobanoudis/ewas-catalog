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
    print("[AUGMENT] Performing comprehensive join (retaining all CpGs)...")
    df_mappings_clean = df_mappings[['cpg', 'chr', 'Start_hg38', 'End_hg38', 'gene']].rename(columns={
        'chr': 'cpg_chr',
        'Start_hg38': 'cpg_start',
        'End_hg38': 'cpg_end'
    })
    
    df_living['input'] = df_living['input'].astype(str).str.strip()
    df_augmented = pd.merge(df_living, df_mappings_clean, left_on='input', right_on='gene', how='right')
    
    if 'gene' in df_augmented.columns:
        df_augmented = df_augmented.drop(columns=['gene'])

    # 2. EWAS Atlas Integration & Gap Analysis
    if Path(atlas_path).is_file():
        print("[AUGMENT] Processing EWAS Atlas data and gap analysis...")
        df_atlas = pd.read_csv(atlas_path, low_memory=False)
        
        # A. Attach traits (with refined formatting/scoring)
        df_augmented = attach_ewas_atlas_traits(df_augmented, df_atlas)
        
        # B. Update ewas_atlas.csv on disk with rank_score
        def calc_rank_score(row):
            try:
                # Only calculate if rank is NOT blank
                if pd.notna(row['rank']) and pd.notna(row['total_associations']) and row['total_associations'] != 0:
                    score = float(row['rank']) / float(row['total_associations'])
                    return f"{score:.3f}"
            except:
                pass
            return "" # Return empty string for CSV if no rank exists
        
        df_atlas_disk = df_atlas.copy()
        df_atlas_disk['rank_score'] = df_atlas_disk.apply(calc_rank_score, axis=1)
        
        # Clean up .0 decimals for CSV readability
        for col in ['rank', 'total_associations', 'pmid']:
            if col in df_atlas_disk.columns:
                # Convert to numeric, then to string without .0, preserving NaN
                df_atlas_disk[col] = pd.to_numeric(df_atlas_disk[col], errors='coerce')
                df_atlas_disk[col] = df_atlas_disk[col].apply(lambda x: f"{int(x)}" if pd.notna(x) else "")


        # Reorder columns to put rank_score between p and rank
        cols = list(df_atlas_disk.columns)
        if 'p' in cols and 'rank' in cols:
            cols.remove('rank_score')
            p_idx = cols.index('p')
            cols.insert(p_idx + 1, 'rank_score')
        df_atlas_disk = df_atlas_disk[cols]
        df_atlas_disk.to_csv(atlas_path, index=False)
        print(f"[AUGMENT] Updated {atlas_path} with rank_score column and integer strings.")

        # C. Analyze unmapped genes (passing original unexpanded mappings)
        unmapped_genes, unmapped_regions, unaccounted_df, appendix_c_df = analyze_unmapped_genes(df_atlas, df_living, df_mappings)
        
        # D. Attach columns
        df_augmented['ewas_unmapped_gene'] = df_augmented['cpg'].map(unmapped_genes)
        df_augmented['ewas_unmapped_regions'] = df_augmented['cpg'].map(unmapped_regions)
        
        # E. Export diagnostic tables
        print("[AUGMENT] Exporting gap analysis diagnostic tables...")
        unaccounted_df.to_csv("unaccounted_genes.csv", index=False)
        appendix_c_df.to_csv("appendix_c_source.csv", index=False)
    else:
        print(f"[AUGMENT] Warning: {atlas_path} not found.")

    # 3. Rename existing 'disgenet_evidence' to 'disgenet_psych_evidence'
    if 'disgenet_evidence' in df_augmented.columns:
        print("[AUGMENT] Renaming 'disgenet_evidence' to 'disgenet_psych_evidence'...")
        df_augmented = df_augmented.rename(columns={'disgenet_evidence': 'disgenet_psych_evidence'})

    # 4. Attach Broad-Spectrum DisGeNET
    if Path(gea_path).is_file():
        print("[AUGMENT] Attaching Broad-Spectrum DisGeNET associations...")
        df_gea = pd.read_csv(gea_path)
        
        # We target 'symbol' column
        df_augmented = annotate_with_disgenet(df_augmented, df_gea, psych_only=False)
        
        # Clean up the duplicate broad column (it adds 'disgenet_diseases' while 'disgenet_disease' might exist)
        # Fix for duplicate disgenet columns: ensure only one broad column exists
        if 'disgenet_diseases' in df_augmented.columns and 'disgenet_disease' in df_augmented.columns:
             print("[AUGMENT] Removing duplicate disgenet broad column...")
             df_augmented = df_augmented.drop(columns=['disgenet_diseases'])
        
        # Drop broad evidence if it appeared
        if 'disgenet_evidence' in df_augmented.columns:
             df_augmented = df_augmented.drop(columns=['disgenet_evidence'])
    else:
        print(f"[AUGMENT] Warning: {gea_path} not found.")

    # Reorder columns
    cols = list(df_augmented.columns)
    cpg_cols = ['cpg', 'cpg_chr', 'cpg_start', 'cpg_end']
    for c in reversed(cpg_cols):
        if c in cols:
            cols.remove(c)
            cols.insert(0, c)
    
    # Put unmapped cols near end
    for c in ['ewas_unmapped_gene', 'ewas_unmapped_regions']:
        if c in cols:
            cols.remove(c)
            cols.append(c)
        
    df_augmented = df_augmented[cols]

    print(f"[AUGMENT] Saving augmented file to: {output_path}")
    df_augmented.to_excel(output_path, index=False)

    # 5. Generate Visualizations
    print("[AUGMENT] Generating visualizations...")
    viz_dir = Path("viz")
    viz_dir.mkdir(exist_ok=True)
    
    if 'disgenet_psych_diseases' in df_augmented.columns:
        plot_disease_heatmap(df_augmented, 
                             disease_col='disgenet_psych_diseases', 
                             output_path="viz/heatmap_psychiatric_seaborn.png",
                             title="Psychiatric Gene-Disease Associations (Seaborn)",
                             dpi=300)
    
    # Use disgenet_disease (the correct one) for the network graph
    target_broad_col = 'disgenet_disease' if 'disgenet_disease' in df_augmented.columns else 'disgenet_diseases'
    if target_broad_col in df_augmented.columns:
        plot_broad_spectrum_network(df_augmented, 
                                    disease_col=target_broad_col, 
                                    output_path="viz/network_broad_spectrum.png",
                                    title="Broad-Spectrum Gene-Disease Network",
                                    score_threshold=0.1,
                                    wrap_width=20)

    if shutil.which("Rscript"):
        r_script_path = Path("disgenet/plot_heatmap.R")
        if r_script_path.is_file():
            print("[AUGMENT] Executing R script for DisGeNET2R heatmap...")
            try:
                subprocess.run(["Rscript", str(r_script_path)], check=True)
            except subprocess.CalledProcessError as e:
                print(f"[AUGMENT] Error running R script: {e}")

    # 6. Generate Summary Report
    print("[AUGMENT] Generating summary report...")
    generate_summary_report(df_augmented, output_path="summary_report.md")

    print("[AUGMENT] Done.")

if __name__ == "__main__":
    augment_living_file()

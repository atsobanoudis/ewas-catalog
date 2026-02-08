import pandas as pd
from pathlib import Path
import sys
import os
import subprocess
import shutil

# Add root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from original_annotation.modules.cpg_integration import load_pi_cpg_mappings, attach_ewas_atlas_traits
from original_annotation.modules.disgenet_utils import annotate_with_disgenet
from original_annotation.modules.viz_engine import plot_disease_heatmap, plot_broad_spectrum_network


def augment_living_file(living_path: str = "annotated_genes.xlsx", 
                        mapping_path: str = "ewas_res_groupsig_128.xlsx",
                        atlas_path: str = "data/ewas_atlas.csv",
                        gea_path: str = "data/disgenet_gea.csv",
                        output_path: str = "annotated_genes_augmented.xlsx"):
    """
    Loads the existing living file and appends CpG coordinates, Atlas traits, and broad DisGeNET data.
    """
    # ... (rest of function)

    print(f"[AUGMENT] Loading living file: {living_path}")
    df_living = pd.read_excel(living_path)
    
    print(f"[AUGMENT] Loading CpG mappings from PI file: {mapping_path}")
    df_mappings = load_pi_cpg_mappings(mapping_path)
    
    # 1. Join CpG coordinates to the living file
    # We join on 'symbol' from living file and 'gene' from mappings
    # Note: One gene can have multiple CpGs, so the row count might increase
    print("[AUGMENT] Joining CpG coordinates...")
    # Clean mappings to just what we need for join
    df_mappings_clean = df_mappings[['cpg', 'chr', 'Start_hg38', 'End_hg38', 'gene']].rename(columns={
        'chr': 'cpg_chr',
        'Start_hg38': 'cpg_start',
        'End_hg38': 'cpg_end'
    })
    
    # Use 'input' as the join key as requested, ensuring no leading/trailing whitespace
    df_living['input'] = df_living['input'].astype(str).str.strip()
    df_augmented = pd.merge(df_living, df_mappings_clean, left_on='input', right_on='gene', how='left')
    
    # Drop the extra 'gene' column from join
    if 'gene' in df_augmented.columns:
        df_augmented = df_augmented.drop(columns=['gene'])

    # 2. Attach EWAS Atlas Traits
    if Path(atlas_path).is_file():
        print("[AUGMENT] Attaching EWAS Atlas traits...")
        df_atlas = pd.read_csv(atlas_path, low_memory=False)
        df_augmented = attach_ewas_atlas_traits(df_augmented, df_atlas)
    else:
        print(f"[AUGMENT] Warning: {atlas_path} not found.")

    # Rename existing 'disgenet_evidence' to 'disgenet_psych_evidence' if it exists
    if 'disgenet_evidence' in df_augmented.columns:
        print("[AUGMENT] Renaming 'disgenet_evidence' to 'disgenet_psych_evidence'...")
        df_augmented = df_augmented.rename(columns={'disgenet_evidence': 'disgenet_psych_evidence'})

    # 3. Attach Broad-Spectrum DisGeNET
    if Path(gea_path).is_file():
        print("[AUGMENT] Attaching Broad-Spectrum DisGeNET associations...")
        df_gea = pd.read_csv(gea_path)
        # We use annotate_with_disgenet which expects 'symbol' column
        # Since we might already have disgenet_psych columns, this will add disgenet_diseases (broad)
        df_augmented = annotate_with_disgenet(df_augmented, df_gea, psych_only=False)
        
        # We only want the disease list, not the evidence for broad spectrum
        if 'disgenet_evidence' in df_augmented.columns:
             df_augmented = df_augmented.drop(columns=['disgenet_evidence'])
    else:
        print(f"[AUGMENT] Warning: {gea_path} not found.")

    # Reorder columns to put CpG info first
    cols = list(df_augmented.columns)
    cpg_cols = ['cpg', 'cpg_chr', 'cpg_start', 'cpg_end']
    for c in reversed(cpg_cols):
        if c in cols:
            cols.remove(c)
            cols.insert(0, c)
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

    print("[AUGMENT] Done.")



if __name__ == "__main__":
    augment_living_file()

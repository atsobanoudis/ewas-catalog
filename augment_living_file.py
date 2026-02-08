import pandas as pd
from pathlib import Path
import sys
import os

# Add root to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from original_annotation.modules.cpg_integration import load_pi_cpg_mappings, attach_ewas_atlas_traits
from original_annotation.modules.disgenet_utils import annotate_with_disgenet

def augment_living_file(living_path: str = "annotated_genes.xlsx", 
                        mapping_path: str = "ewas_res_groupsig_128.xlsx",
                        atlas_path: str = "data/ewas_atlas.csv",
                        gea_path: str = "data/disgenet_gea.csv",
                        output_path: str = "annotated_genes_augmented.xlsx"):
    """
    Loads the existing living file and appends CpG coordinates, Atlas traits, and broad DisGeNET data.
    """
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
    
    df_augmented = pd.merge(df_living, df_mappings_clean, left_on='symbol', right_on='gene', how='left')
    
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

    # 3. Attach Broad-Spectrum DisGeNET
    if Path(gea_path).is_file():
        print("[AUGMENT] Attaching Broad-Spectrum DisGeNET associations...")
        df_gea = pd.read_csv(gea_path)
        # We use annotate_with_disgenet which expects 'symbol' column
        # Since we might already have disgenet_psych columns, this will add disgenet_diseases (broad)
        df_augmented = annotate_with_disgenet(df_augmented, df_gea, psych_only=False)
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
    print("[AUGMENT] Done.")

if __name__ == "__main__":
    augment_living_file()

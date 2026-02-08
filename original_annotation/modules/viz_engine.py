import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

def plot_disease_heatmap(df: pd.DataFrame, 
                         disease_col: str = 'disgenet_diseases', 
                         gene_col: str = 'symbol',
                         output_path: str = "disease_heatmap.png",
                         title: str = "Gene-Disease Association Heatmap",
                         cmap: str = "viridis",
                         figsize: tuple = (12, 10),
                         font_scale: float = 1.0):
    """
    Generates a publication-quality heatmap showing associations between genes and diseases.
    Expects associations in the format "Disease Name, Score;\nNext Disease, Score".
    """
    # ... (parsing logic)
    parsed_data = []
    
    for _, row in df.iterrows():
        gene = row[gene_col]
        assoc_str = row[disease_col]
        
        if pd.isna(assoc_str) or assoc_str == 'n/a':
            continue
            
        # Split by semicolon and optional newline
        parts = re.split(r';\n|;', assoc_str)
        for part in parts:
            part = part.strip()
            if not part:
                continue
                
            # Expecting "Disease, Score"
            if ',' in part:
                try:
                    disease, score_str = part.rsplit(',', 1)
                    score = float(score_str.strip())
                    parsed_data.append({
                        'gene': gene,
                        'disease': disease.strip(),
                        'score': score
                    })
                except (ValueError, TypeError):
                    # Skip 'error' or non-numeric scores for heatmap
                    continue

    if not parsed_data:
        print("[VIZ] No numeric association data found for heatmap.")
        return

    viz_df = pd.DataFrame(parsed_data)
    
    # 2. Pivot to matrix form
    pivot_df = viz_df.pivot_table(index='gene', columns='disease', values='score', fill_value=0)
    
    # 3. Plotting (Academic & Clean style)
    sns.set_theme(style="white", font_scale=font_scale)
    plt.figure(figsize=figsize)
    
    # Use specified cmap (default viridis is color-blind friendly)
    ax = sns.heatmap(pivot_df, 
                     annot=False, 
                     cmap=cmap, 
                     linewidths=.5, 
                     cbar_kws={'label': 'DisGeNET Association Score'})
    
    plt.title(title, fontsize=16 * font_scale, pad=20)
    plt.xlabel("Disease", fontsize=12 * font_scale)
    plt.ylabel("Gene", fontsize=12 * font_scale)
    
    # Rotate labels for readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    # Save output
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"[VIZ] Heatmap saved to: {output_path}")
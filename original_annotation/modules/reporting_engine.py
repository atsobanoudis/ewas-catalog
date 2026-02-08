import pandas as pd
import re
from datetime import datetime

def generate_summary_report(df: pd.DataFrame, output_path: str = "summary_report.md"):
    """
    Generates a summary report highlighting top gene-disease associations.
    """
    
    # 1. Analyze Top Psychiatric Hits
    psych_hits = []
    if 'disgenet_psych_diseases' in df.columns:
        for _, row in df.iterrows():
            gene = row.get('symbol', 'Unknown')
            assoc_str = row['disgenet_psych_diseases']
            
            if pd.isna(assoc_str) or assoc_str == 'n/a':
                continue
                
            parts = re.split(r';\n|;', assoc_str)
            for part in parts:
                if ',' in part:
                    try:
                        disease, score_str = part.rsplit(',', 1)
                        score = float(score_str.strip())
                        if score >= 0.3: # Filter for significance in report
                            psych_hits.append({
                                'Gene': gene,
                                'Disease': disease.strip(),
                                'Score': score
                            })
                    except (ValueError, TypeError):
                        continue
    
    # Sort by score descending
    psych_hits.sort(key=lambda x: x['Score'], reverse=True)
    top_psych_hits = psych_hits[:10] # Top 10

    # 2. Analyze Top EWAS Atlas Hits (based on frequency?)
    # Since Atlas doesn't have a score in the same way, we can list genes with most traits?
    atlas_hits = []
    if 'ewas_atlas_traits' in df.columns:
        for _, row in df.iterrows():
            gene = row.get('symbol', 'Unknown')
            traits_str = row['ewas_atlas_traits']
            if pd.isna(traits_str) or traits_str == 'n/a':
                continue
            
            traits = [t.strip() for t in re.split(r';\n|;', traits_str) if t.strip()]
            if traits:
                atlas_hits.append({
                    'Gene': gene,
                    'Trait_Count': len(traits),
                    'Traits': ", ".join(traits[:3]) + ("..." if len(traits) > 3 else "")
                })
    
    atlas_hits.sort(key=lambda x: x['Trait_Count'], reverse=True)
    top_atlas_hits = atlas_hits[:10]

    # 3. Generate Markdown Content
    lines = []
    lines.append(f"# Psychiatric Epigenetics Discovery Report")
    lines.append(f"**Date:** {datetime.now().strftime('%Y-%m-%d')}\n")
    
    lines.append("## Executive Summary")
    lines.append("This report synthesizes key findings from the annotated gene list, highlighting significant psychiatric disease associations and epigenetic trait overlaps.\n")
    
    lines.append("## Top Psychiatric Associations (DisGeNET)")
    lines.append("Genes with the strongest evidence for psychiatric disorders (Score >= 0.3):")
    lines.append("| Gene | Disease | Score |")
    lines.append("|---|---|---|")
    if top_psych_hits:
        for hit in top_psych_hits:
            lines.append(f"| {hit['Gene']} | {hit['Disease']} | {hit['Score']} |")
    else:
        lines.append("| None | No significant hits found | - |")
    lines.append("\n")
    
    lines.append("## Top Epigenetic Trait Clusters (EWAS Atlas)")
    lines.append("Genes with the highest number of associated traits in EWAS studies:")
    lines.append("| Gene | Trait Count | Top Traits |")
    lines.append("|---|---|---|")
    if top_atlas_hits:
        for hit in top_atlas_hits:
            lines.append(f"| {hit['Gene']} | {hit['Trait_Count']} | {hit['Traits']} |")
    else:
        lines.append("| None | No Atlas data found | - |")
    lines.append("\n")
    
    lines.append("## Interpretation")
    lines.append(f"Identified {len(psych_hits)} significant psychiatric gene-disease associations.")
    if top_psych_hits:
        top_gene = top_psych_hits[0]['Gene']
        top_disease = top_psych_hits[0]['Disease']
        lines.append(f"The strongest signal is observed for **{top_gene}** in association with **{top_disease}**.")
    
    # Write to file
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    
    print(f"[REPORT] Summary report written to: {output_path}")
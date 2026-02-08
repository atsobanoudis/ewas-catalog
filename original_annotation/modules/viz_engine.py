import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
import networkx as nx
import textwrap

def plot_disease_heatmap(df: pd.DataFrame, 
                         disease_col: str = 'disgenet_diseases', 
                         gene_col: str = 'symbol',
                         output_path: str = "disease_heatmap.png",
                         title: str = "Gene-Disease Association Heatmap",
                         cmap: str = "viridis",
                         figsize: tuple = None,
                         font_scale: float = 1.0,
                         dpi: int = 300):
    """
    Generates a publication-quality heatmap showing associations between genes and diseases.
    Expects associations in the format "Disease Name, Score;\nNext Disease, Score".
    """
    # 1. Parse the associations into a long-form DataFrame
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
    
    # Dynamically adjust figsize if not provided
    if figsize is None:
        width = max(12, pivot_df.shape[1] * 0.3)
        height = max(10, pivot_df.shape[0] * 0.3)
        figsize = (width, height)

    plt.figure(figsize=figsize)
    
    # Create a mask for 0 values to make them appear blank
    mask = pivot_df == 0
    
    # Use specified cmap (default viridis is color-blind friendly)
    ax = sns.heatmap(pivot_df, 
                     mask=mask,
                     annot=False, 
                     cmap=cmap, 
                     linewidths=.5, 
                     cbar_kws={'label': 'DisGeNET Association Score'},
                     facecolor='white') # Blank cells are white
    
    plt.title(title, fontsize=16 * font_scale, pad=20)
    plt.xlabel("Disease", fontsize=12 * font_scale)
    plt.ylabel("Gene", fontsize=12 * font_scale)
    
    # Rotate labels for readability
    plt.xticks(rotation=45, ha='right', fontsize=8 * font_scale)
    plt.yticks(rotation=0, fontsize=8 * font_scale)
    
    plt.tight_layout()
    
    # Save output with high DPI for zooming
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    print(f"[VIZ] Heatmap saved to: {output_path} (DPI: {dpi})")

def plot_broad_spectrum_network(df: pd.DataFrame, 
                                disease_col: str = 'disgenet_diseases', 
                                gene_col: str = 'symbol',
                                output_path: str = "broad_spectrum_network.png",
                                title: str = "Gene-Disease Network",
                                score_threshold: float = 0.1,
                                wrap_width: int = 15):
    """
    Generates a network graph for broad-spectrum gene-disease associations.
    Genes are source nodes, Diseases are target nodes.
    Filters out associations with score < score_threshold.
    Wraps disease labels.
    """
    G = nx.Graph()
    
    for _, row in df.iterrows():
        gene = row[gene_col]
        assoc_str = row[disease_col]
        
        if pd.isna(assoc_str) or assoc_str == 'n/a':
            continue
            
        # We will add the gene node later only if it has valid edges
        
        parts = re.split(r';\n|;', assoc_str)
        for part in parts:
            part = part.strip()
            if not part:
                continue
            
            if ',' in part:
                try:
                    disease, score_str = part.rsplit(',', 1)
                    disease = disease.strip()
                    score = float(score_str.strip())
                    
                    if score >= score_threshold:
                        # Wrap disease name
                        disease_label = "\n".join(textwrap.wrap(disease, wrap_width))
                        
                        # Add nodes and edge
                        G.add_node(gene, type='gene', label=gene)
                        G.add_node(disease, type='disease', label=disease_label)
                        G.add_edge(gene, disease, weight=score)
                        
                except (ValueError, TypeError):
                    continue

    if G.number_of_nodes() == 0:
        print(f"[VIZ] No nodes to plot in network graph (threshold={score_threshold}).")
        return

    plt.figure(figsize=(20, 20)) # Larger figure
    
    # Layout - Kamada Kawai often handles disjoint components nicely
    try:
        pos = nx.kamada_kawai_layout(G)
    except:
        pos = nx.spring_layout(G, k=0.3, iterations=50) # Fallback
    
    # Nodes
    genes = [n for n, d in G.nodes(data=True) if d.get('type') == 'gene']
    diseases = [n for n, d in G.nodes(data=True) if d.get('type') == 'disease']
    
    # Genes (Blue)
    nx.draw_networkx_nodes(G, pos, nodelist=genes, node_color='skyblue', node_size=300, alpha=0.9, label='Genes')
    
    # Diseases (Red) - Sized by degree? Or fixed? Let's keep fixed but distinct color
    nx.draw_networkx_nodes(G, pos, nodelist=diseases, node_color='salmon', node_size=150, alpha=0.7, label='Diseases')
    
    # Edges - Width based on weight
    edges = G.edges(data=True)
    weights = [d['weight'] for u, v, d in edges]
    # Normalize weights for width (e.g., 0.1 to 1.0 -> width 0.5 to 5)
    widths = [w * 5 for w in weights]
    
    nx.draw_networkx_edges(G, pos, edge_color='gray', width=widths, alpha=0.4)
    
    # Labels
    # Use the 'label' attribute we set (wrapped for diseases)
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_color='black')
    
    plt.title(f"{title} (Score > {score_threshold})", fontsize=24)
    plt.axis('off')
    
    # Legend
    plt.legend(scatterpoints=1)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"[VIZ] Network graph saved to: {output_path}")
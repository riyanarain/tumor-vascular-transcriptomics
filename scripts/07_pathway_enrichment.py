"""
Script 06: Pathway Enrichment Analysis (GSEA)
Project: Tumor Vascular Transcriptomics
Author: Riya Narain
Date: March 2026

Goal: Identify which biological pathways are enriched in Cluster 2 (Tumor Angiogenic ECs)
      compared to Cluster 1 (Normal Quiescent ECs) using pre-ranked GSEA.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
import warnings
warnings.filterwarnings('ignore')

# Install gseapy if needed:
# pip install gseapy

import gseapy as gp
from gseapy import barplot, dotplot

# ─────────────────────────────────────────────
# 1. Load DEG results
# ─────────────────────────────────────────────
print("Loading DEG results...")

deg = pd.read_csv('../results/cluster2_vs_cluster1_DEG.csv')
print(f"  Loaded {len(deg)} genes")
print(deg.head())

# ─────────────────────────────────────────────
# 2. Prepare ranked gene list for prerank GSEA
#    Ranking metric: sign(logFC) * -log10(pval)
#    This captures both direction and significance
# ─────────────────────────────────────────────
print("\nPreparing ranked gene list...")

deg = deg.dropna(subset=['logfoldchanges', 'pvals'])
deg['pvals_safe'] = deg['pvals'].clip(lower=1e-300)  # avoid log(0)
deg['rank_metric'] = np.sign(deg['logfoldchanges']) * -np.log10(deg['pvals_safe'])
deg = deg.sort_values('rank_metric', ascending=False)
deg = deg.drop_duplicates(subset='gene', keep='first')

# Filter out ribosomal protein genes (RPL/RPS) before ranking
ribo_mask = ~deg['gene'].str.match(r'^(RPL|RPS|MRPL|MRPS)')
deg = deg[ribo_mask]
print(f"  After removing ribosomal genes: {len(deg)} genes")


ranked_genes = deg[['gene', 'rank_metric']].set_index('gene')
print(f"  Ranked {len(ranked_genes)} unique genes")
print(f"  Top upregulated: {ranked_genes.head(5).index.tolist()}")
print(f"  Top downregulated: {ranked_genes.tail(5).index.tolist()}")

# Save ranked list
os.makedirs('results', exist_ok=True)
ranked_genes.to_csv('../results/ranked_gene_list.csv')
print("  Saved: results/ranked_gene_list.csv")

# ─────────────────────────────────────────────
# 3. Run prerank GSEA — KEGG pathways
# ─────────────────────────────────────────────
print("\nRunning prerank GSEA on KEGG pathways...")

os.makedirs('../results/gsea_kegg', exist_ok=True)

kegg_results = gp.prerank(
    rnk=ranked_genes,
    gene_sets='KEGG_2021_Human',
    outdir='results/gsea_kegg',
    min_size=10,
    max_size=500,
    permutation_num=1000,
    seed=42,
    verbose=False
)

kegg_df = kegg_results.res2d
kegg_df = kegg_df.sort_values('NES', ascending=False)
kegg_df.to_csv('../results/gsea_kegg_results.csv', index=False)
print(f"  KEGG: {len(kegg_df)} pathways tested")
print(f"  Significant (FDR<0.25): {(kegg_df['FDR q-val'] < 0.25).sum()}")

# ─────────────────────────────────────────────
# 4. Run prerank GSEA — GO Biological Process
# ─────────────────────────────────────────────
print("\nRunning prerank GSEA on GO Biological Process...")

os.makedirs('results/gsea_gobp', exist_ok=True)

gobp_results = gp.prerank(
    rnk=ranked_genes,
    gene_sets='GO_Biological_Process_2023',
    outdir='../results/gsea_gobp',
    min_size=10,
    max_size=500,
    permutation_num=1000,
    seed=42,
    verbose=False
)

gobp_df = gobp_results.res2d
gobp_df = gobp_df.sort_values('NES', ascending=False)
gobp_df.to_csv('../results/gsea_gobp_results.csv', index=False)
print(f"  GO BP: {len(gobp_df)} pathways tested")
print(f"  Significant (FDR<0.25): {(gobp_df['FDR q-val'] < 0.25).sum()}")

# ─────────────────────────────────────────────
# 5. Visualize top results — KEGG
# ─────────────────────────────────────────────
print("\nGenerating figures...")

os.makedirs('../figures/pathway_enrichment', exist_ok=True)

def plot_top_pathways(df, title, filename, n=15, fdr_cutoff=0.25):
    """Bar plot of top enriched and depleted pathways."""
    sig = df[df['FDR q-val'] < fdr_cutoff].copy()
    if len(sig) == 0:
        print(f"  Warning: No significant pathways at FDR<{fdr_cutoff} for {title}")
        sig = df.head(20)  # show top 20 by NES anyway

    top_up = sig[sig['NES'] > 0].head(n)
    top_down = sig[sig['NES'] < 0].tail(n)
    plot_df = pd.concat([top_up, top_down]).sort_values('NES')

    # Shorten long pathway names
    plot_df = plot_df.copy()
    plot_df['Term'] = plot_df['Term'].str.replace(r'\(GO:\d+\)', '', regex=True).str.strip()
    plot_df['Term'] = plot_df['Term'].apply(lambda x: x[:60] + '...' if len(x) > 60 else x)

    colors = ['#d73027' if x > 0 else '#4575b4' for x in plot_df['NES']]

    fig, ax = plt.subplots(figsize=(10, max(6, len(plot_df) * 0.4)))
    bars = ax.barh(range(len(plot_df)), plot_df['NES'], color=colors, edgecolor='white', linewidth=0.5)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['Term'], fontsize=9)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)

    # Add FDR labels
    for i, (nes, fdr) in enumerate(zip(plot_df['NES'], plot_df['FDR q-val'])):
        label = f'FDR={fdr:.3f}'
        x_pos = nes + (0.05 if nes > 0 else -0.05)
        ha = 'left' if nes > 0 else 'right'
        ax.text(x_pos, i, label, va='center', ha=ha, fontsize=7, color='#333333')

    legend_elements = [
        plt.Rectangle((0,0),1,1, color='#d73027', label='Enriched in Tumor ECs (Cluster 2)'),
        plt.Rectangle((0,0),1,1, color='#4575b4', label='Enriched in Normal ECs (Cluster 1)')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

plot_top_pathways(
    kegg_df,
    'KEGG Pathway Enrichment: Tumor vs Normal Endothelial Cells',
    '../figures/pathway_enrichment/kegg_enrichment_barplot.png'
)

plot_top_pathways(
    gobp_df,
    'GO Biological Process Enrichment: Tumor vs Normal Endothelial Cells',
    '../figures/pathway_enrichment/gobp_enrichment_barplot.png'
)

# ─────────────────────────────────────────────
# 6. Enrichment plot for key pathways of interest
# ─────────────────────────────────────────────
def plot_enrichment_curve(results_obj, term_keyword, filename):
    """Plot GSEA enrichment curve for a specific pathway."""
    terms = results_obj.res2d['Term'].tolist()
    match = [t for t in terms if term_keyword.lower() in t.lower()]
    if not match:
        print(f"  Warning: No pathway matching '{term_keyword}' found")
        return
    term = match[0]
    try:
        ax = results_obj.plot(terms=term, show_ranking=True, legend_kws={'loc': 'upper left'})
        plt.title(f'GSEA Enrichment: {term[:60]}', fontsize=10, fontweight='bold')
        plt.tight_layout()
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {filename}")
    except Exception as e:
        print(f"  Could not plot enrichment curve for '{term}': {e}")

plot_enrichment_curve(kegg_results, 'angiogen', '../figures/pathway_enrichment/enrichment_curve_angiogenesis.png')
plot_enrichment_curve(kegg_results, 'ECM', '../figures/pathway_enrichment/enrichment_curve_ecm.png')
plot_enrichment_curve(gobp_results, 'blood vessel', '../figures/pathway_enrichment/enrichment_curve_blood_vessel.png')

# ─────────────────────────────────────────────
# 7. Summary table of top findings
# ─────────────────────────────────────────────
print("\nTop KEGG pathways enriched in Tumor ECs (Cluster 2):")
top_kegg = kegg_df[kegg_df['NES'] > 0].head(10)[['Term', 'NES', 'FDR q-val', 'Lead_genes']]
print(top_kegg.to_string(index=False))

print("\nTop GO BP pathways enriched in Tumor ECs (Cluster 2):")
top_gobp = gobp_df[gobp_df['NES'] > 0].head(10)[['Term', 'NES', 'FDR q-val', 'Lead_genes']]
print(top_gobp.to_string(index=False))

print("\n✅ Pathway enrichment analysis complete!")
print("Outputs saved to:")
print("  results/gsea_kegg_results.csv")
print("  results/gsea_gobp_results.csv")
print("  results/ranked_gene_list.csv")
print("  figures/pathway_enrichment/")

"""
Remove immune contamination (Cluster 6), add biological annotations to clusters,
and save the clean final dataset for downstream analysis and sharing.
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import os
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 1

# Load clustered endothelial dataset
print("Loading clustered endothelial dataset...")
adata = sc.read_h5ad('../data/processed/endothelial_clustered.h5ad')

# Remove Cluster 5
contamination_mask = adata.obs['leiden'] != '6'
adata_clean = adata[contamination_mask].copy()

# Add biological cluster annotations, based on differential expression results and marker gene analysis
cluster_annotations = {
    '0': 'Venular/Stalk EC',        # ACKR1, CLU, C7, VWF
    '1': 'Normal Quiescent EC',     # FCN3, IL7R, EPAS1
    '2': 'Tumor Angiogenic EC',     # RGCC, SPARC, COL4A1, CD34, PLVAP 
    '3': 'Lymphatic EC',            # CCL21, TFF3, PROX1, PDPN
    '4': 'Normal Activated EC',     # HPGD, ITM2B, EDNRB
    '5': 'Arterial EC',             # GJA5, EFNB2, IGFBP3
}

adata_clean.obs['cell_state'] = adata_clean.obs['leiden'].map(cluster_annotations)
adata_clean.obs['cell_state'] = adata_clean.obs['cell_state'].astype('category')

print("  Cluster annotations:")
for k, v in cluster_annotations.items():
    n = (adata_clean.obs['leiden'] == k).sum()
    pct_tumor = (
        adata_clean.obs.loc[adata_clean.obs['leiden'] == k, 'tissue_category'] == 'Tumor'
    ).mean() * 100
    print(f"    Cluster {k} → {v:30s} (n={n}, {pct_tumor:.0f}% tumor)")

# Add a tumor-specific flag for Cluster 2
adata_clean.obs['is_tumor_angiogenic'] = (
    adata_clean.obs['leiden'] == '2'
).astype(bool)

# Save final annotated dataset
output_path = '../data/processed/endothelial_final_annotated.h5ad'
adata_clean.write_h5ad(output_path)

# Generate final summary figures
os.makedirs('../figures/final', exist_ok=True)

# Color palette for biological clusters
cluster_colors = {
    'Venular/Stalk EC':      '#e41a1c',
    'Normal Quiescent EC':   '#377eb8',
    'Tumor Angiogenic EC':   '#ff7f00',
    'Lymphatic EC':          '#4daf4a',
    'Normal Activated EC':   '#984ea3',
    'Arterial EC':           '#a65628',
}

# Create UMAP colored by biological cell state
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

sc.pl.umap(
    adata_clean,
    color='cell_state',
    palette=cluster_colors,
    title='Endothelial Cell Subtypes\n(Lung Adenocarcinoma)',
    legend_loc='right margin',
    ax=axes[0],
    show=False,
    frameon=False
)

sc.pl.umap(
    adata_clean,
    color='tissue_category',
    palette={'Tumor': '#d73027', 'Normal': '#4575b4'},
    title='Tissue Origin\n(Tumor vs Normal)',
    legend_loc='right margin',
    ax=axes[1],
    show=False,
    frameon=False
)

plt.suptitle('Final Annotated Endothelial Cell Atlas — LUAD', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('../figures/final/umap_final_annotated.png', dpi=150, bbox_inches='tight')
plt.close()

# Make composition bar chart with percentage tumor per cluster
fig, ax = plt.subplots(figsize=(9, 5))

comp = adata_clean.obs.groupby(['cell_state', 'tissue_category']).size().unstack(fill_value=0)
comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100
comp_pct = comp_pct.reindex(list(cluster_annotations.values()))

comp_pct.plot(
    kind='bar',
    stacked=True,
    color=['#d73027', '#4575b4'],
    ax=ax,
    edgecolor='white',
    linewidth=0.5
)

ax.set_xlabel('')
ax.set_ylabel('Percentage of Cells (%)', fontsize=12)
ax.set_title('Tumor vs Normal Composition per EC Subtype', fontsize=13, fontweight='bold')
ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right', fontsize=10)
ax.legend(title='Tissue', fontsize=10)
ax.axhline(50, color='black', linestyle='--', linewidth=0.8, alpha=0.5, label='50%')

# annotate Cluster 2 bar
for i, label in enumerate(comp_pct.index):
    if label == 'Tumor Angiogenic EC':
        ax.annotate('★ Key Finding', xy=(i, 120), ha='center', fontsize=9,
                    color='#ff7f00', fontweight='bold')

plt.tight_layout()
plt.savefig('../figures/final/composition_barplot.png', dpi=150, bbox_inches='tight')
plt.close()

# Export final cell metadata as csv
metadata = adata_clean.obs[['leiden', 'cell_state', 'tissue_category', 'is_tumor_angiogenic']].copy()
metadata.to_csv('../results/final_cell_metadata.csv')
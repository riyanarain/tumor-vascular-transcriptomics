import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("../data/processed/endothelial_clustered.h5ad")

# Differential expression between clusters
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# Get top marker genes for each cluster
sc.pl.rank_genes_groups(adata, n_genes=20)

# Save results
result = adata.uns['rank_genes_groups']
genes_df = pd.DataFrame({
    group: result['names'][group][:50] 
    for group in result['names'].dtype.names
})
genes_df.to_csv("../results/cluster_marker_genes.csv")

"""
Visualize differential expression results
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path

# Setup
figures_dir = Path("../figures/DE_analysis")
figures_dir.mkdir(parents=True, exist_ok=True)

# Load data
adata = sc.read_h5ad("../data/processed/endothelial_clustered.h5ad")

# Remove Cluster 6 (immune contamination)
adata_clean = adata[adata.obs['leiden'] != '6'].copy()

print("="*60)
print("DIFFERENTIAL EXPRESSION VISUALIZATION")
print("="*60)

# IMPORTANT: Rerun differential expression on cleaned data
print("\nRunning differential expression on cleaned data...")
sc.tl.rank_genes_groups(adata_clean, groupby='leiden', method='wilcoxon')
print("✓ DE analysis complete")

# Figure 1: Dotplot of top markers per cluster
print("\n1. Creating dotplot of top markers...")

top_markers = {
    '0': ['ACKR1', 'VWF', 'CLU', 'C7'],
    '1': ['FCN3', 'EPAS1', 'IL7R', 'CAV1'],
    '2': ['RGCC', 'SPARC', 'COL4A1', 'CD34', 'PLVAP'],
    '3': ['CCL21', 'TFF3', 'PROX1', 'PDPN'],
    '4': ['HPGD', 'EDNRB', 'HLA-E', 'ITM2B'],
    '5': ['GJA5', 'EFNB2', 'IGFBP3', 'CXCL12']
}

# Flatten to list
all_markers = []
for cluster, genes in top_markers.items():
    all_markers.extend(genes)

# Filter for genes that exist
available_markers = [g for g in all_markers if g in adata_clean.var_names]

fig, ax = plt.subplots(figsize=(12, 6))
sc.pl.dotplot(adata_clean, available_markers, groupby='leiden', 
              ax=ax, show=False, dendrogram=False)
plt.title('Top Marker Genes per Cluster', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / 'dotplot_top_markers.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Dotplot saved")

# Figure 2: Heatmap of top markers
print("\n2. Creating heatmap...")

fig = sc.pl.rank_genes_groups_heatmap(
    adata_clean, 
    n_genes=10, 
    groupby='leiden',
    show_gene_labels=True,
    show=False,
    save='_marker_heatmap.png'
    
)
print("✓ Heatmap saved")

# Figure 3: Violin plots for key Cluster 2 markers
print("\n3. Creating violin plots for Cluster 2 markers...")

cluster2_markers = ['RGCC', 'SPARC', 'PLVAP', 'COL4A1', 'CD34']
available_cluster2 = [m for m in cluster2_markers if m in adata_clean.var_names]

fig, axes = plt.subplots(1, len(available_cluster2), figsize=(4*len(available_cluster2), 4))
if len(available_cluster2) == 1:
    axes = [axes]

for idx, gene in enumerate(available_cluster2):
    sc.pl.violin(adata_clean, gene, groupby='leiden', ax=axes[idx], show=False)
    axes[idx].set_title(f'{gene}', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / 'violin_cluster2_markers.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Violin plots saved")

# Figure 4: UMAP showing key markers
print("\n4. Creating UMAP feature plots...")

feature_genes = ['PLVAP', 'SPARC', 'GJA5', 'CCL21', 'ACKR1']
available_features = [g for g in feature_genes if g in adata_clean.var_names]

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, gene in enumerate(available_features):
    sc.pl.umap(adata_clean, color=gene, ax=axes[idx], show=False, 
               cmap='Reds', vmax='p99', frameon=False)
    axes[idx].set_title(f'{gene} Expression', fontsize=12, fontweight='bold')

# Hide unused subplots
for idx in range(len(available_features), 6):
    axes[idx].axis('off')

plt.tight_layout()
plt.savefig(figures_dir / 'umap_marker_expression.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ UMAP feature plots saved")

# Figure 5: Stacked violin - Tumor vs Normal for Cluster 2 markers
print("\n5. Creating tumor vs normal comparison...")

fig, axes = plt.subplots(1, len(available_cluster2), figsize=(4*len(available_cluster2), 4))
if len(available_cluster2) == 1:
    axes = [axes]

for idx, gene in enumerate(available_cluster2):
    sc.pl.violin(adata_clean, gene, groupby='tissue_category', ax=axes[idx], show=False)
    axes[idx].set_title(f'{gene}: Tumor vs Normal', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(figures_dir / 'violin_tumor_vs_normal.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Tumor vs normal comparison saved")

# Figure 6: Focused comparison - Cluster 2 vs Cluster 1
print("\n6. Running focused DE: Cluster 2 vs Cluster 1...")

sc.tl.rank_genes_groups(
    adata_clean,
    groupby='leiden',
    groups=['2'],
    reference='1',
    method='wilcoxon',
    key_added='de_cluster2_vs_cluster1'
)

# Extract results
result = adata_clean.uns['de_cluster2_vs_cluster1']
cluster2_vs_1_genes = pd.DataFrame({
    'gene': result['names']['2'],
    'logfoldchanges': result['logfoldchanges']['2'],
    'pvals': result['pvals']['2'],
    'pvals_adj': result['pvals_adj']['2']
})

# Filter significant
sig_genes = cluster2_vs_1_genes[cluster2_vs_1_genes['pvals_adj'] < 0.05].copy()
print(f"\nSignificant genes (Cluster 2 vs 1): {len(sig_genes)}")

# Save
sig_genes.to_csv("../results/cluster2_vs_cluster1_DEG.csv", index=False)

# Show top upregulated in Cluster 2
print("\nTop 20 upregulated in Cluster 2 (Tumor Angiogenic):")
top20 = sig_genes.nlargest(20, 'logfoldchanges')[['gene', 'logfoldchanges', 'pvals_adj']]
print(top20.to_string(index=False))

# Volcano plot
fig, ax = plt.subplots(figsize=(10, 8))

# All genes
ax.scatter(cluster2_vs_1_genes['logfoldchanges'], 
          -np.log10(cluster2_vs_1_genes['pvals_adj'] + 1e-300),  # Add small value to avoid log(0)
          c='gray', alpha=0.5, s=10)

# Significant upregulated
sig_up = sig_genes[sig_genes['logfoldchanges'] > 0]
ax.scatter(sig_up['logfoldchanges'], 
          -np.log10(sig_up['pvals_adj'] + 1e-300),
          c='#e74c3c', alpha=0.7, s=20, label='Upregulated in Cluster 2')

# Significant downregulated
sig_down = sig_genes[sig_genes['logfoldchanges'] < 0]
ax.scatter(sig_down['logfoldchanges'], 
          -np.log10(sig_down['pvals_adj'] + 1e-300),
          c='#3498db', alpha=0.7, s=20, label='Downregulated in Cluster 2')

# Label top genes
top_genes_to_label = sig_genes.nlargest(10, 'logfoldchanges')['gene'].tolist()
for _, row in sig_genes.iterrows():
    if row['gene'] in top_genes_to_label:
        ax.text(row['logfoldchanges'], -np.log10(row['pvals_adj'] + 1e-300), 
               row['gene'], fontsize=8, alpha=0.8)

ax.axhline(-np.log10(0.05), color='black', linestyle='--', linewidth=1, alpha=0.5)
ax.set_xlabel('Log2 Fold Change', fontsize=12)
ax.set_ylabel('-log10(Adjusted P-value)', fontsize=12)
ax.set_title('Volcano Plot: Cluster 2 (Tumor Angiogenic) vs Cluster 1 (Normal Quiescent)', 
            fontsize=14, fontweight='bold')
ax.legend()
plt.tight_layout()
plt.savefig(figures_dir / 'volcano_cluster2_vs_cluster1.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Volcano plot saved")

print("\n" + "="*60)
print("ALL FIGURES SAVED!")
print("="*60)
print(f"\nOutput directory: {figures_dir}/")
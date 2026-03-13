"""
Exploratory analysis of endothelial cells from GSE131907.
Generates PCA, UMAP, highly variable genes plot, and angiogenesis marker dotplot.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd

# Setup
figures_dir = Path("../figures/meeting")
figures_dir.mkdir(parents=True, exist_ok=True)

# Load in endothelial data
adata = sc.read_h5ad("../data/processed/endothelial_cells.h5ad")
tissue_counts = adata.obs['tissue_category'].value_counts()

# Make cell distribution bar chart
fig, ax = plt.subplots(figsize=(8, 6))
tissue_counts.plot(kind='bar', ax=ax, color=['#e74c3c', '#3498db'])
ax.set_title('Endothelial Cell Distribution: Tumor vs Normal', fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue Type', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(figures_dir / 'cell_distribution.png', dpi=300, bbox_inches='tight')
plt.close()

# Make figure for sample breakdown
fig, ax = plt.subplots(figsize=(10, 6))
sample_counts = adata.obs.groupby(['Sample_Origin', 'tissue_category']).size().unstack(fill_value=0)
sample_counts.plot(kind='bar', ax=ax, stacked=True, color=['#e74c3c', '#3498db'])
ax.set_title('Endothelial Cells by Sample Origin', fontsize=14, fontweight='bold')
ax.set_xlabel('Sample Origin', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.legend(title='Tissue Category')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(figures_dir / 'sample_breakdown.png', dpi=300, bbox_inches='tight')
plt.close()

# Run PCA for dimensionality reduction 
sc.tl.pca(adata)
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.pca(adata, color='tissue_category', ax=ax, show=False, 
          palette={'Tumor': '#e74c3c', 'Normal': '#3498db'})
plt.title('PCA: Tumor vs Normal Endothelial Cells', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / 'pca_analysis.png', dpi=300, bbox_inches='tight')
plt.close()

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
n_variable = adata.var['highly_variable'].sum()
sc.pl.highly_variable_genes(adata, show=False)
plt.suptitle(f'Highly Variable Genes (n={n_variable})', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / 'variable_genes.png', dpi=300, bbox_inches='tight')
plt.close()

# Run UMAP
sc.pp.neighbors(adata, n_neighbors=15)  # run with default n_neighbors = 15 
sc.tl.umap(adata)
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(adata, color='tissue_category', ax=ax, show=False,
           palette={'Tumor': '#e74c3c', 'Normal': '#3498db'})
plt.title('UMAP: Tumor vs Normal Endothelial Cells', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / 'umap_analysis.png', dpi=300, bbox_inches='tight')
plt.close()

# Check presence of known angiogenesis markers 
markers = ['VEGFA', 'PECAM1', 'VWF', 'CDH5', 'FLT1', 'KDR', 'ANGPT2']  # from literature 
available_markers = [m for m in markers if m in adata.var_names]

if len(available_markers) >= 3:
    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.dotplot(adata, available_markers, groupby='tissue_category', ax=ax, show=False)
    plt.title('Angiogenesis Marker Expression', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(figures_dir / 'marker_expression.png', dpi=300, bbox_inches='tight')
    plt.close()

# Save summary statistics
summary = {
    'Total Cells': adata.n_obs,
    'Total Genes': adata.n_vars,
    'Tumor ECs': (adata.obs['tissue_category'] == 'Tumor').sum(),
    'Normal ECs': (adata.obs['tissue_category'] == 'Normal').sum(),
    'Highly Variable Genes': adata.var['highly_variable'].sum() if 'highly_variable' in adata.var else 'N/A',
    'Available Markers': ', '.join(available_markers)
}

summary_df = pd.DataFrame(summary.items(), columns=['Metric', 'Value'])
summary_df.to_csv(figures_dir / 'summary_statistics.csv', index=False)
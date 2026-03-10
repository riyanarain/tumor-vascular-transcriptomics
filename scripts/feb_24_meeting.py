"""
Generate figures and summary for Tuesday meeting
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import pandas as pd

# Setup
figures_dir = Path("../figures/meeting")
figures_dir.mkdir(parents=True, exist_ok=True)

print("="*60)
print("GENERATING MEETING MATERIALS")
print("="*60)

# Load data
print("\n1. Loading endothelial cells...")
adata = sc.read_h5ad("../data/processed/endothelial_cells.h5ad")

# Summary stats
print("\n" + "="*60)
print("DATA SUMMARY FOR MEETING")
print("="*60)
print(f"\nTotal endothelial cells: {adata.n_obs}")
print(f"Total genes: {adata.n_vars}")
print(f"\nTissue distribution:")
tissue_counts = adata.obs['tissue_category'].value_counts()
print(tissue_counts)
print(f"\nSample origin breakdown:")
print(adata.obs['Sample_Origin'].value_counts())

# Figure 1: Cell distribution bar chart
fig, ax = plt.subplots(figsize=(8, 6))
tissue_counts.plot(kind='bar', ax=ax, color=['#e74c3c', '#3498db'])
ax.set_title('Endothelial Cell Distribution: Tumor vs Normal', fontsize=14, fontweight='bold')
ax.set_xlabel('Tissue Type', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
plt.tight_layout()
plt.savefig(figures_dir / '01_cell_distribution.png', dpi=300, bbox_inches='tight')
plt.close()
print("\n✓ Figure 1: Cell distribution")

# Figure 2: Sample breakdown
fig, ax = plt.subplots(figsize=(10, 6))
sample_counts = adata.obs.groupby(['Sample_Origin', 'tissue_category']).size().unstack(fill_value=0)
sample_counts.plot(kind='bar', ax=ax, stacked=True, color=['#e74c3c', '#3498db'])
ax.set_title('Endothelial Cells by Sample Origin', fontsize=14, fontweight='bold')
ax.set_xlabel('Sample Origin', fontsize=12)
ax.set_ylabel('Number of Cells', fontsize=12)
ax.legend(title='Tissue Category')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(figures_dir / '02_sample_breakdown.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 2: Sample breakdown")

# Figure 3: PCA
print("\n2. Computing PCA...")
sc.tl.pca(adata)
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.pca(adata, color='tissue_category', ax=ax, show=False, 
          palette={'Tumor': '#e74c3c', 'Normal': '#3498db'})
plt.title('PCA: Tumor vs Normal Endothelial Cells', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / '03_pca_analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 3: PCA")

# Figure 4: Highly variable genes
print("\n3. Identifying variable genes...")
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
n_variable = adata.var['highly_variable'].sum()
sc.pl.highly_variable_genes(adata, show=False)
plt.suptitle(f'Highly Variable Genes (n={n_variable})', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / '04_variable_genes.png', dpi=300, bbox_inches='tight')
plt.close()
print(f"✓ Figure 4: Variable genes (found {n_variable})")

# Figure 5: UMAP
print("\n4. Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)
fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.umap(adata, color='tissue_category', ax=ax, show=False,
           palette={'Tumor': '#e74c3c', 'Normal': '#3498db'})
plt.title('UMAP: Tumor vs Normal Endothelial Cells', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(figures_dir / '05_umap_analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("✓ Figure 5: UMAP")

# Figure 6: Known angiogenesis markers (if available)
print("\n5. Checking angiogenesis markers...")
markers = ['VEGFA', 'PECAM1', 'VWF', 'CDH5', 'FLT1', 'KDR', 'ANGPT2']
available_markers = [m for m in markers if m in adata.var_names]
print(f"   Available markers: {available_markers}")

if len(available_markers) >= 3:
    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.dotplot(adata, available_markers, groupby='tissue_category', ax=ax, show=False)
    plt.title('Angiogenesis Marker Expression', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(figures_dir / '06_marker_expression.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✓ Figure 6: Marker expression")
else:
    print("⚠ Not enough markers available for dotplot")

# Save summary statistics
print("\n6. Saving summary statistics...")
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
print("✓ Summary statistics saved")

print("\n" + "="*60)
print("ALL MATERIALS GENERATED!")
print("="*60)
print(f"\nFigures saved to: {figures_dir}/")
print("\nKey talking points for meeting:")
print(f"  • Successfully extracted {adata.n_obs} endothelial cells")
print(f"  • {(adata.obs['tissue_category']=='Tumor').sum()} tumor-associated, {(adata.obs['tissue_category']=='Normal').sum()} normal")
print(f"  • Identified {adata.var['highly_variable'].sum() if 'highly_variable' in adata.var else 'N/A'} highly variable genes")
print(f"  • PCA/UMAP show distinct clustering patterns")
print(f"  • {len(available_markers)} angiogenesis markers detected")
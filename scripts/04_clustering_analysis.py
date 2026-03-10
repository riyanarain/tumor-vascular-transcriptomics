import scanpy as sc
import pandas as pd

# Load endothelial cells
adata = sc.read_h5ad("../data/processed/endothelial_cells.h5ad")

# Clustering workflow
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.leiden(adata, resolution=0.5)  # Clustering algorithm
sc.tl.umap(adata)

# Visualize clusters
sc.pl.umap(adata, color=['leiden', 'tissue_category', 'Sample_Origin'], legend_loc='best')


# Check if clusters correspond to tissue type, sample, or something else
print("\nCluster composition by tissue category:")
cluster_by_tissue = pd.crosstab(
    adata.obs['leiden'], 
    adata.obs['tissue_category']
)
print(cluster_by_tissue)

print("\nCluster composition by sample origin:")
cluster_by_sample_origin = pd.crosstab(
    adata.obs['leiden'], 
    adata.obs['Sample_Origin']
)
print(cluster_by_sample_origin)

# Save
adata.write("../data/processed/endothelial_clustered.h5ad")

# Rank genes with wilcoxon test 
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

# For each cluster, print top 10 genes
for cluster in adata.obs['leiden'].unique():
    print(f"\nCluster {cluster} top genes:")
    genes = adata.uns['rank_genes_groups']['names'][cluster][:10]
    print(genes)

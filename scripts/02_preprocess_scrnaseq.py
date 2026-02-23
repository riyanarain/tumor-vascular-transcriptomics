"""
Preprocess single-cell RNA-seq data (Memory-efficient version)
Load, filter, and normalize GSE131907 data
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import anndata

# Set up paths
data_dir = Path("../data/raw/GSE131907")
output_dir = Path("../data/processed")
output_dir.mkdir(exist_ok=True)

print("Loading annotation data...")
annotations = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    sep="\t",
    index_col="Index"
)

print(f"Loaded annotations for {len(annotations)} cells")
print(f"\nCell types found:")
print(annotations['Cell_type'].value_counts())

print(f"\nSample origins:")
print(annotations['Sample_Origin'].value_counts())

# Use raw UMI matrix instead - it's smaller and we can normalize ourselves
print("\nLoading raw UMI matrix (more memory-efficient)...")
print("This may take 2-3 minutes...")

adata = sc.read_text(
    data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    delimiter='\t',
    first_column_names=True
).T  # Transpose so cells are rows

print(f"\nLoaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Add annotations - match cell IDs
common_cells = adata.obs_names.intersection(annotations.index)
print(f"Matching cells between expression and annotation: {len(common_cells)}")

# Subset to matching cells only
adata = adata[common_cells, :]
adata.obs = annotations.loc[common_cells]

print(f"\nAfter matching: {adata.n_obs} cells × {adata.n_vars} genes")

# Basic QC filtering
print("\nApplying QC filters...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")

# Normalize and log-transform
print("\nNormalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("\nSaving processed data...")
adata.write(output_dir / "GSE131907_processed.h5ad")

print("\n✓ Preprocessing complete!")
print(f"Processed data saved to: {output_dir / 'GSE131907_processed.h5ad'}")
print(f"\nData summary:")
print(f"  Cells: {adata.n_obs}")
print(f"  Genes: {adata.n_vars}")
print(f"  Endothelial cells: {(adata.obs['Cell_type'] == 'Endothelial cells').sum()}")
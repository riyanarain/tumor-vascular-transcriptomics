"""
Preprocess single-cell RNA-seq data
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

# Using raw UMI matrix because smaller 
adata = sc.read_text(
    data_dir / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    delimiter='\t',
    first_column_names=True
).T  # transpose so cells are rows

print(f"\nLoaded: {adata.n_obs} cells × {adata.n_vars} genes")

# Add annotations
common_cells = adata.obs_names.intersection(annotations.index)

# Subset to matching cells only
adata = adata[common_cells, :]
adata.obs = annotations.loc[common_cells]

# QC filtering
print("\nApplying QC filters...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write(output_dir / "GSE131907_processed.h5ad")
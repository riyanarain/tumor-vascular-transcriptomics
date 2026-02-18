"""
Preprocess single-cell RNA-seq data
Load, filter, and normalize GSE131907 data
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

# Set up paths
data_dir = Path("../data/raw/GSE131907")
output_dir = Path("../data/processed")
output_dir.mkdir(exist_ok=True)

print("Loading annotation data...")
# Load cell annotations
annotations = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    sep="\t"
)

print(f"Loaded annotations for {len(annotations)} cells")
print(f"\nCell types found:")
print(annotations['Cell_type'].value_counts())

print(f"\nSample origins:")
print(annotations['Sample_Origin'].value_counts())

print("\nLoading expression matrix (this may take a few minutes)...")
# Load expression matrix - genes as rows, cells as columns
expr_matrix = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
    sep="\t",
    index_col=0
)

print(f"Expression matrix: {expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} cells")

# Create AnnData object (transpose so cells are rows, genes are columns)
print("\nCreating AnnData object...")
adata = sc.AnnData(X=expr_matrix.T)

# Add cell annotations
# Match by Index column
annotations_indexed = annotations.set_index('Index')
adata.obs = annotations_indexed.loc[adata.obs_names]

print(f"\nFinal AnnData object: {adata.n_obs} cells × {adata.n_vars} genes")

# Basic filtering
print("\nApplying basic QC filters...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells × {adata.n_vars} genes")

# Save processed data
print("\nSaving processed data...")
adata.write(output_dir / "GSE131907_processed.h5ad")

print("\n✓ Preprocessing complete!")
print(f"Processed data saved to: {output_dir / 'GSE131907_processed.h5ad'}")
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

print("Loading data...")

# Load expression matrix (you'll need to download this manually)
# Option 1: If you have the .txt.gz file
adata = sc.read_text(data_dir / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz")

# Load cell annotations
annotations = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    sep="\t"
)

# Add annotations to adata
adata.obs = adata.obs.join(annotations.set_index('cell_id'))

print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes")

# Basic filtering (adjust thresholds as needed)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(f"After filtering: {adata.n_obs} cells and {adata.n_vars} genes")

# Save processed data
adata.write(output_dir / "GSE131907_processed.h5ad")

print("Preprocessing complete!")
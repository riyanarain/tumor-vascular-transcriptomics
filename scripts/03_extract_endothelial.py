"""
Extract endothelial cells from tumor vs normal tissue
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# Load processed data
adata = sc.read_h5ad("../data/processed/GSE131907_processed.h5ad")

print(f"Total cells: {adata.n_obs}")
print("\nCell type distribution:")
print(adata.obs['cell_type'].value_counts())

# Extract endothelial cells
# (Check the actual column names in the annotation file)
endothelial = adata[adata.obs['cell_type'].str.contains('Endothelial', case=False, na=False)]

print(f"\nExtracted {endothelial.n_obs} endothelial cells")

# Separate tumor vs normal
tumor_ec = endothelial[endothelial.obs['tissue_type'] == 'Tumor']
normal_ec = endothelial[endothelial.obs['tissue_type'] == 'Normal']

print(f"Tumor ECs: {tumor_ec.n_obs}")
print(f"Normal ECs: {normal_ec.n_obs}")

# Save subsets
endothelial.write("../data/processed/endothelial_cells.h5ad")
tumor_ec.write("../data/processed/tumor_endothelial.h5ad")
normal_ec.write("../data/processed/normal_endothelial.h5ad")

print("Endothelial cell extraction complete!")
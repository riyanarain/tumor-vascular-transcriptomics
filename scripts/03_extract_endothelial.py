"""
Extract endothelial cells from tumor vs normal tissue
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# Load in processed data
adata = sc.read_h5ad("../data/processed/GSE131907_processed.h5ad")

# Extract just endothelial cells
endothelial = adata[adata.obs['Cell_type'] == 'Endothelial cells'].copy()
print(f"Extracted {endothelial.n_obs} endothelial cells")

# Categorize tissue types
tumor_samples = ['tLung', 'tL/B', 'mBrain', 'mLN', 'PE']  # tumor-related samples
normal_samples = ['nLung', 'nLN']  # normal tissue

endothelial.obs['tissue_category'] = endothelial.obs['Sample_Origin'].apply(
    lambda x: 'Tumor' if x in tumor_samples else 'Normal'
)

# Separate tumor vs normal for easy access
tumor_ec = endothelial[endothelial.obs['tissue_category'] == 'Tumor'].copy()
normal_ec = endothelial[endothelial.obs['tissue_category'] == 'Normal'].copy()

print(f"\nTumor-associated ECs: {tumor_ec.n_obs}")
print(f"Normal ECs: {normal_ec.n_obs}")

# Save all subsets 
output_dir = Path("../data/processed")
endothelial.write(output_dir / "endothelial_cells.h5ad", compression="gzip")
tumor_ec.write(output_dir / "tumor_endothelial.h5ad", compression="gzip")
normal_ec.write(output_dir / "normal_endothelial.h5ad", compression="gzip")
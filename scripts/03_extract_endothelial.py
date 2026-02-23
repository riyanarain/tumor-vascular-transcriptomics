"""
Extract endothelial cells from tumor vs normal tissue
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

# Load processed data
print("Loading processed data...")
adata = sc.read_h5ad("../data/processed/GSE131907_processed.h5ad")

print(f"Total cells: {adata.n_obs}")
print("\nCell type distribution:")
print(adata.obs['Cell_type'].value_counts())

# Extract endothelial cells
print("\n" + "="*50)
print("Extracting endothelial cells...")
endothelial = adata[adata.obs['Cell_type'] == 'Endothelial cells'].copy()

print(f"Extracted {endothelial.n_obs} endothelial cells")

# Categorize tissue types
print("\nCategorizing tissue origin...")
tumor_samples = ['tLung', 'tL/B', 'mBrain', 'mLN', 'PE']  # tumor-related samples
normal_samples = ['nLung', 'nLN']  # normal tissue

endothelial.obs['tissue_category'] = endothelial.obs['Sample_Origin'].apply(
    lambda x: 'Tumor' if x in tumor_samples else 'Normal'
)

print("\nSample origin breakdown:")
print(endothelial.obs['Sample_Origin'].value_counts())

print("\nTissue category:")
print(endothelial.obs['tissue_category'].value_counts())

# Separate tumor vs normal for easy access
tumor_ec = endothelial[endothelial.obs['tissue_category'] == 'Tumor'].copy()
normal_ec = endothelial[endothelial.obs['tissue_category'] == 'Normal'].copy()

print(f"\nTumor-associated ECs: {tumor_ec.n_obs}")
print(f"Normal ECs: {normal_ec.n_obs}")

# Save all versions
output_dir = Path("../data/processed")
print("\nSaving endothelial cell data...")

endothelial.write(output_dir / "endothelial_cells.h5ad", compression="gzip")
tumor_ec.write(output_dir / "tumor_endothelial.h5ad", compression="gzip")
normal_ec.write(output_dir / "normal_endothelial.h5ad", compression="gzip")

print("\n✓ Endothelial cell extraction complete!")
print(f"\nFiles saved:")
print(f"  - All ECs: {output_dir / 'endothelial_cells.h5ad'}")
print(f"  - Tumor ECs: {output_dir / 'tumor_endothelial.h5ad'}")
print(f"  - Normal ECs: {output_dir / 'normal_endothelial.h5ad'}")
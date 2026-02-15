"""
Download GSE131907 data from GEO
Author: Noemi Reche-Ley, Riya Narain
Date: 2026-02-14
"""

import GEOparse
import os
from pathlib import Path

# Setup directories
data_dir = Path("../data/raw/GSE131907")
data_dir.mkdir(parents=True, exist_ok=True)

print("Downloading GSE131907 from GEO...")

# Download the dataset
gse = GEOparse.get_GEO(geo="GSE131907", destdir=str(data_dir))

print(f"Dataset: {gse.name}")
print(f"Title: {gse.metadata['title'][0]}")
print(f"Number of samples: {len(gse.gsms)}")

# Download supplementary files
print("\nDownloading supplementary files...")
# Note: You may need to manually download large files from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907

print("""
Manual download needed for large files:
1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907
2. Download:
   - GSE131907_Lung_Cancer_cell_annotation.txt.gz
   - GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz
   - GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
3. Place in: data/raw/GSE131907/
""")

print("\nData download script complete!")
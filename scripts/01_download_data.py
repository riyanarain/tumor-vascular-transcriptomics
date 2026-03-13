"""
Download GSE131907 data from GEO
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

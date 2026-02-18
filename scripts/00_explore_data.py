"""
Explore the structure of GSE131907 data files
"""

import pandas as pd
from pathlib import Path

# Use the correct path where files actually are
data_dir = Path("../data/raw/GSE131907")

print("=== Checking Annotation File ===")
annotations = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    sep="\t",
    nrows=5  # Just first 5 rows to peek
)

print("\nColumn names:")
print(annotations.columns.tolist())

print("\nFirst few rows:")
print(annotations.head())

print("\n" + "="*50)
print("=== Checking Expression Matrix ===")

# Peek at expression matrix structure
expr = pd.read_csv(
    data_dir / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz",
    sep="\t",
    nrows=5,
    index_col=0
)

print("\nExpression matrix shape (first 5 rows):")
print(f"Genes: {expr.shape[0]}, Cells: {expr.shape[1]}")

print("\nFirst few cell IDs:")
print(expr.columns[:10].tolist())

print("\nFirst few gene names:")
print(expr.index[:10].tolist())
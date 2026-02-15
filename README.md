# Tumor Vascular Transcriptomics

Integrated transcriptomic and morphometric analysis of tumor-induced vascular dysfunction in lung adenocarcinoma.

## Project Overview
This project analyzes single-cell RNA-seq data (GSE131907) to identify molecular signatures of tumor-associated endothelial dysfunction and integrates these findings with morphometric analysis of vascularized microphysiological tumor platforms (µPTM).

## Team
- Noemi Reche-Ley
- Riya Narain

## Datasets
- **GSE131907**: Single-cell RNA-seq from 44 lung adenocarcinoma patients (208,506 cells)
- **µPTM Images**: Confocal microscopy of perfusable microvascular networks

## Workflow
1. Download and preprocess GSE131907 data
2. Extract and annotate endothelial cells
3. Perform differential expression analysis (tumor vs. normal ECs)
4. Pathway enrichment analysis (GSEA)
5. Morphometric feature extraction from images
6. Machine learning integration of transcriptomic and morphometric data

## Requirements
- R (≥4.0) with Seurat, DESeq2, clusterProfiler
- Python (≥3.8) with pandas, scikit-learn, scanpy
- ImageJ/Fiji with Angiogenesis Analyzer plugin

## Installation
```bash
# Clone repository
git clone https://github.com/YOUR-USERNAME/tumor-vascular-transcriptomics.git

# Install R packages
Rscript -e "install.packages(c('Seurat', 'DESeq2', 'clusterProfiler'))"

# Install Python packages
pip install pandas numpy scikit-learn scanpy matplotlib seaborn
```

## Usage
[To be added as scripts are developed]

## Results
[To be added]

## References
- Kim N et al. (2020) Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. Nat Commun.
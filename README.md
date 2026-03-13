# Tumor Vascular Transcriptomics
Integrated transcriptomic analysis of tumor-induced vascular dysfunction in lung adenocarcinoma.

**BMDS 205 · Stanford University · Winter 2026**
**Authors:** Noemi Reche-Ley & Riya Narain

---

## Project Overview
We applied unsupervised clustering and differential expression analysis to single-cell RNA sequencing data from 44 lung adenocarcinoma patients to identify transcriptionally distinct endothelial cell populations. We identified a tumor-enriched angiogenic EC cluster (Cluster 2) defined by PLVAP, CD34, COL4A1, SPARC, and RGCC, consistent with increased vascular permeability and ECM remodeling.

**Dataset:** GSE131907 (Kim et al., Nat Commun 2020) — 208,506 cells from 44 LUAD patients

---

## Repository Structure
```
scripts/          # Analysis pipeline (run in order)
figures/          # All generated figures
results/          # CSV outputs from analyses
data/processed/   # Small processed files (large files excluded)
```

---

## Analysis Pipeline

| Script | Description |
|--------|-------------|
| `01_download_data.py` | Download GSE131907 from GEO |
| `02_preprocess_scrnaseq.py` | QC, CPM normalization, log transform |
| `03_extract_endothelial.py` | Extract 1,996 endothelial cells |
| `04_exploratory_analysis.py` | PCA, UMAP, marker validation |
| `05_clustering_analysis.py` | Leiden clustering, cluster annotation |
| `06_differential_expression.py` | Wilcoxon DE, volcano plot |
| `07_pathway_enrichment.py` | GSEA against KEGG and GO BP |
| `08_final_cleanup.py` | Remove contamination, annotate final dataset |

---

## Installation
```bash
# Clone repository
git clone https://github.com/riyanarain/tumor-vascular-transcriptomics.git
cd tumor-vascular-transcriptomics

# Create conda environment
conda create -n tumor-vascular python=3.10
conda activate tumor-vascular

# Install dependencies
pip install -r requirements.txt
```

---

## Data
Raw data files are not included due to size. Download from GEO:
- **GSE131907**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907
- Run `scripts/01_download_data.py` to download automatically

---

## Key Findings
- 6 biologically distinct EC populations identified
- Cluster 2 (Tumor Angiogenic EC): 59% tumor-derived, 1.8-fold enrichment
- Key markers: PLVAP (+2.35 log₂FC), CD34 (+2.22), COL4A1 (+1.78), SPARC, RGCC
- GSEA confirms ECM remodeling and biosynthetic reprogramming in tumor ECs
- Normal ECs enriched for immune surveillance pathways

---

## Requirements
See `requirements.txt` for full list. Key packages:
- scanpy >= 1.9.0
- gseapy >= 1.0.4
- lifelines >= 0.27.0
- pandas, numpy, matplotlib, seaborn

---

## Reference
Kim N et al. Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma. *Nat Commun* 11, 2285 (2020).

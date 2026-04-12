# Immune Infiltration and Survival in Breast Cancer (TCGA BRCA)

This project aims to investigate whether the level of immune cell infiltration in breast cancer tumors is associated with patient survival outcomes, using publicly available data from The Cancer Genome Atlas (TCGA).

---

## Background

The tumor immune microenvironment plays a critical role in cancer progression and patient outcomes. Tumors that attract more immune cells — particularly cytotoxic T cells and NK cells — are generally associated with better prognosis, while immunosuppressive cell types like regulatory T cells (Tregs) can dampen anti-tumor immunity.

This analysis offers to shed light on if breast cancer patients with higher immune cell infiltration live longer**

To answer this, immune infiltration scores were calculated for five cell types using mRNA expression data, and Kaplan-Meier survival analysis was used to compare outcomes between patients with high vs. low infiltration.

---

## Data

- **Source: (https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018) — TCGA Breast Invasive Carcinoma (BRCA) dataset
- **Expression data:** mRNA sequencing (RSEM normalized), 1,082 patients
- **Clinical data:** Overall survival (OS_MONTHS), disease-specific survival status (DSS_STATUS), 1,084 patients
- **Final merged cohort:** 1,062 patients after quality filtering

---

## Methods

### Survival endpoint
A composite cancer-death variable (`tumor_death`) was derived by combining overall survival status (OS_STATUS) and disease-specific survival status (DSS_STATUS), to distinguish confirmed cancer-related deaths from non-cancer deaths and censored patients.

### Immune scoring
Immune infiltration scores were calculated for five cell types using curated gene signatures. For each cell type, expression values were log2-transformed (log2(x+1)) and averaged across the signature genes per patient:

| Cell Type | Genes |
|-----------|-------|
| CD8 T cells | CD8A, CD8B, GZMB, PRF1, IFNG |
| Tregs | FOXP3, IL2RA, CTLA4 |
| B cells | CD19, MS4A1, CD79A |
| NK cells | NCAM1, NKG7, KLRD1 |
| Macrophages | CD68, CD163, MRC1 |

### Survival analysis
Patients were stratified into high and low infiltration groups using the median score for each cell type. Kaplan-Meier survival curves were plotted for each group and compared using the log-rank test.

---

## Key Findings

### Kaplan-Meier Survival Curves

![Kaplan-Meier curves](km_curves.png)

Three cell types showed statistically significant associations with survival:

| Cell Type | p-value | Direction |
|-----------|---------|-----------|
| CD8 T cells | 0.007 | High = better survival |
| B cells | 0.004 | High = better survival |
| NK cells | 0.047 | High = better survival |
| Tregs | 0.495 | Not significant |
| Macrophages | 0.715 | Not significant |

CD8 T cells and B cells showed the strongest survival associations, consistent with their roles as key mediators of adaptive anti-tumor immunity. NK cells showed a borderline significant effect, reflecting their cytotoxic function in innate immunity.

Tregs showed no significant survival difference despite being immunosuppressive — a finding that likely reflects the known tendency of immune-inflamed tumors to co-recruit both effector and regulatory populations. Macrophages seemed to show no survival signal, consistent with their functional heterogeneity in the tumor microenvironment (pro- vs. anti-tumor subtypes).

### Immune Score Correlation Heatmap

![Correlation heatmap](immune_correlation_heatmap.png)

All five immune scores were positively correlated, with the strongest correlation observed between CD8 T cells and Tregs (r=0.85). This counterintuitive finding is consistent with published literature — tumors that are broadly immune-inflamed tend to recruit all immune populations simultaneously, rather than selectively attracting either effector or suppressive cells.

Macrophages showed the weakest correlations overall, particularly with B cells (r=0.46), reflecting their distinct recruitment mechanisms as innate immune cells.

---

## Tools & Libraries

- Python 3
- pandas, numpy
- lifelines (Kaplan-Meier, log-rank test)
- matplotlib, seaborn

---

## Repository Structure

```
brca-survival-analysis/
├── data/                        # Raw data files (not included, download from cBioPortal)
├── proj.py                      # Main analysis script
├── km_curves.png                # Kaplan-Meier survival curves
├── immune_correlation_heatmap.png  # Immune score correlation heatmap
└── README.md
```

---

## How to Run

1. Download TCGA BRCA data from [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_tcga)
   - `data_mrna_seq_v2_rsem.txt`
   - `data_clinical_patient.txt`
2. Place files in a `data/` folder
3. Update file paths in `proj.py` if needed
4. Install dependencies: `pip install pandas numpy lifelines matplotlib seaborn`
5. Run: `python3 proj.py`


# 📊 Differential Gene Expression Analysis of Human Pancreatic Alpha and Beta Cells

This repository contains all R scripts and resources used to perform a full pipeline for differential gene expression (DGE) analysis using RNA-seq data from the **Human Islets (HPAP)** project. The aim of this project is to identify transcriptional alterations between pancreatic **beta** and **alpha** cells in **type 2 diabetes (T2D)** donors compared to non-diabetic individuals. This analysis is part of the Bachelor Thesis (TFG) by Carla Batlle Simó (UVic-UCC, 2025).

The project integrates public transcriptomic datasets with open-source statistical tools to:
- Process and normalize bulk and pseudobulk RNA-seq data
- Identify differentially expressed genes using DESeq2
- Visualize results (PCA, UMAP, volcano plots, heatmaps)
- Extract meaningful biological insights through enrichment analysis

---

## 📥 Data Acquisition and Format

Data used in this project was downloaded from the [Human Islets portal](https://www.humanislets.com/#/). The data is publicly available but may require prior registration or data request.

### Required files:
- `donor.csv`: clinical metadata per donor
- `unproc_rnaseq.csv`: raw count matrix (bulk RNA-seq)
- `unproc_pbrna_Beta.csv`: raw counts (beta cells pseudobulk)
- `unproc_pbrna_Alpha.csv`: raw counts (alpha cells pseudobulk)

Ensure all files are placed in your R working directory (e.g., `~/Documents/practiques/data`). Metadata must include `record_id`, `donorage`, `diagnosis_computed`, and `donorsex` columns.

---

## 🧬 Project Scripts

### `run_deg_analysis.R`
Runs the DEG analysis between conditions using DESeq2. Includes filtering, normalization, and differential testing. You can adjust filtering parameters inside this script.

### `visualize_deg_results.R`
Visualizes DEG output: volcano plots, heatmaps, PCA and UMAP embeddings, top gene scatterplots.

### `counts2tpm.R`
Function to convert raw gene counts to TPM using gene length annotations.

### `pca_donant_analysis.R`
Performs PCA on donor metadata variables (age, BMI, diagnosis, etc.).

### `full_expression_boxplots.R`
Creates boxplots for global gene expression and individual genes of interest (e.g., INS, GCG, MKI67).

### `distrib_dades.R`
Performs quality control and filters samples and genes based on expression thresholds.

### `get_deg_results.R`
Annotates DEG results with SYMBOL and ENTREZID gene identifiers.

---

## 🚀 Analysis Pipeline

### Step-by-step:
```r
pca_donant_analysis("donor.csv")
full_expression_boxplots("unproc_rnaseq.csv", "TPM_merged.csv")
distrib_bulk <- distrib_dades("unproc_rnaseq.csv")
distrib_beta <- distrib_dades("unproc_pbrna_Beta.csv")
distrib_alpha <- distrib_dades("unproc_pbrna_Alpha.csv")

# DEG Analysis
bulk_results  <- run_deg_analysis("unproc_rnaseq.csv", "donor.csv")
beta_results  <- run_deg_analysis("unproc_pbrna_Beta.csv", "donor.csv")
alpha_results <- run_deg_analysis("unproc_pbrna_Alpha.csv", "donor.csv")

# Visualization
visualize_deg_results(bulk_results,  tag = "BULK")
visualize_deg_results(beta_results,  tag = "BETA")
visualize_deg_results(alpha_results, tag = "ALPHA")

# Annotate results
get_deg_results(beta_results$res_deg)
get_deg_results(alpha_results$res_deg)
```

---

## ⚙️ Parameters and Settings
The following thresholds and options are used by default (modifiable):
- Adjusted p-value (padj) threshold: `0.05`
- Log2 fold change (log2FC) threshold: `1` (absolute)
- Group comparison: `None` vs `Type2`
- Genes used in plots: filtered DEGs only

All thresholds and filtering criteria are adjustable directly in the R scripts:
- Age filters: `donorage >= 10` or `donorage >= 18`
- Zero-expression filter: max 80% of zeros per condition
- Exclude samples manually: e.g., `c("R032", "R036", "R233", ...)`

---

## 🧠 Functional Enrichment (Optional)
You can take the list of DEG symbols and perform enrichment analysis using external tools:

- [Enrichr](https://maayanlab.cloud/Enrichr/): for GO terms, KEGG, WikiPathways
- [WikiPathways](https://www.wikipathways.org/)

Copy the output gene lists (e.g., `clean_list_beta`) and paste them into Enrichr. You may use the full list, or just upregulated/downregulated genes. Outputs can reveal transcription factors (via ChEA/ENCODE), enriched pathways, or disease links.

---

## 🤝 Contact and Contributions
If you find any issues, have questions, or would like to collaborate, feel free to open an issue or submit a pull request. This project is open to improvement and extension.

---

## Author
**Carla Batlle Simó**  
Bachelor's Degree in Biotechnology - UVic-UCC (2025)

This README is part of her Bachelor's Thesis and designed to ensure reproducibility and reusability for other researchers.

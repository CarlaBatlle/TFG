
# Differential Gene Expression Analysis of Human Pancreatic Alpha and Beta Cells

This repository contains all R scripts and resources used to perform a full pipeline for differential gene expression (DGE) analysis using RNA-seq data from the **Human Islets (HPAP)** project. The aim of this project is to identify transcriptional alterations between pancreatic **beta** and **alpha** cells in **type 2 diabetes (T2D)** donors compared to non-diabetic individuals. This analysis is part of the Bachelor Thesis (TFG) by Carla Batlle Simó (UVic-UCC, 2025).

The project integrates public transcriptomic datasets with open-source statistical tools to:
- Process and normalize bulk and pseudobulk RNA-seq data
- Identify differentially expressed genes using DESeq2
- Visualize results (PCA, UMAP, volcano plots, heatmaps)
- Extract meaningful biological insights through enrichment analysis

The workflow includes the next steps:

<img width="275" alt="image" src="https://github.com/user-attachments/assets/f8eb346f-b741-4dd5-a8cc-0f75d9a52207" />

---

## Data Acquisition and Format

Data used in this project was downloaded from the [Human Islets portal](https://www.humanislets.com/#/). The data is publicly available but may require prior registration or data request.

### Download procedure:

1. Go to [Human Islets portal](https://www.humanislets.com/#/)
2. Navigate to the **Data Download** section.
3. Apply the following filter:
   - Diagnosis: **SA (non-diabetic)** and **Type 2 Diabetes**
4. Select the following datasets:
   - **Bulk RNA-seq**
   - **Pseudobulk Alpha**
   - **Pseudobulk Beta**
5. Download the following files:
   - `donor.csv`: clinical metadata per donor
   - `unproc_rnaseq.csv`: raw count matrix (bulk RNA-seq)
   - `unproc_pbrna_Beta.csv`: raw counts matrix (beta cells pseudobulk)
   - `unproc_pbrna_Alpha.csv`: raw counts matrix (alpha cells pseudobulk)

### Required files:
- `donor.csv`: clinical metadata per donor
- `unproc_rnaseq.csv`: raw count matrix (bulk RNA-seq)
- `unproc_pbrna_Beta.csv`: raw counts matrix (beta cells pseudobulk)
- `unproc_pbrna_Alpha.csv`: raw counts matrix (alpha cells pseudobulk)

Ensure all files are placed in your R working directory (e.g., `~/Documents/practiques/data`). 
You also have the data abailable in the repository.

---

##  Project Scripts

### `counts2tpm.R`
Function to convert raw gene counts to TPM using gene length annotations.

### `DEG.R`
Runs the DEG analysis between conditions using DESeq2. Includes data anlysis, filtering, normalization, and differential expressed genes testing. You can adjust filtering parameters inside this script.

---

##  Analysis Pipeline

### Step-by-step:
```r
# Metadata analysis
pca_donant_analysis("donor.csv")

# Outliers filtering
full_expression_boxplots("unproc_rnaseq.csv", "TPM_merged.csv")

# Espression distribution
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

##  Parameters and Settings
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

## Functional Enrichment Analysis

You can take the list of DEG symbols and perform enrichment analysis using external tools:

- [Enrichr](https://maayanlab.cloud/Enrichr/): for GO terms, KEGG, WikiPathways...

  The ones used in the bachelo's thesis:
- [WikiPathways](https://www.wikipathways.org/): for pathway enrichment
- [ChEA/ENCODE (via Enrichr)](https://maayanlab.cloud/Enrichr/#libraries): for transcription factor target enrichment

Copy the output gene lists (e.g., `clean_list_beta`) and paste them into Enrichr.  
You may use the full list, or just upregulated/downregulated genes.  
Outputs can reveal transcription factors (via ChEA/ENCODE), enriched pathways, or disease links.

---

##  Contact and Contributions
If you find any issues, have questions, or would like to collaborate, feel free to open an issue or submit a pull request. This project is open to improvement and extension.

---

## Author
**Carla Batlle Simó**  
Bachelor's Degree in Biotechnology - UVic-UCC (2025)

This README is part of her Bachelor's Thesis and designed to ensure reproducibility and reusability for other researchers.

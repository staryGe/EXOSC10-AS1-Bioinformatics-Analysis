*For Chinese version of this documentation, please see [README_CN.md](./README_CN.md).*

# EXOSC10-AS1 Knockdown Transcriptomic & Functional Enrichment Analysis

![R](https://img.shields.io/badge/R-4.6%2B-blue)
![DESeq2](https://img.shields.io/badge/Bioconductor-DESeq2-green)
![clusterProfiler](https://img.shields.io/badge/Bioconductor-clusterProfiler-orange)

## Project Overview
This project provides a systematic Differential Expression Analysis (DEG) and functional enrichment mining for RNA-seq transcriptomic data (**GSE289213**) targeting lncRNA **EXOSC10-AS1** knockdown in colorectal cancer cell lines (HCT116).

By comparing **NC (Control)** against **siRNA-knockdown** groups, combined with **GO (BP/CC/MF)** and **KEGG Pathway (Pathview)** analyses, this study systematically unveils the downstream molecular mechanisms and pathways modulated by `EXOSC10-AS1` silencing.

---



## Key Biological Insights

1. **Differential Expression Profile (DEGs)**:
   * In the `siRNA_1 vs NC` comparison, **283 significantly altered genes** were identified ($p.adjust < 0.05, \vert{}\log_2 FC\vert{} > 1$), comprising **193 down-regulated** and **90 up-regulated** genes.
2. **Cellular Functional Inactivation (GO Enrichment)**:
   * **BP (Biological Process)**: Primarily suppressed immune cell differentiation and cytokine production.
   * **CC (Cellular Component)**: Markedly disrupted cytoplasmic vesicle lumen structures.
   * **MF (Molecular Function)**: Impaired cytokine receptor binding activities.
3. **Canonical Pathway Response (KEGG Signaling)**:
   * Enriched significantly in the **`Cytokine-cytokine receptor interaction` (hsa04060)** pathway.
   * Pathway mapping via `pathview` confirmed that knockdown specifically down-regulated **TNF receptor family members (e.g., TNFR2, 4-1BB)** and **IL27/IL35** regulatory networks, blocking extracellular signal communication.

---



## Key Results & Visualizations

| Analysis Module            |                        Visual Output                         | Description                                                  |
| :------------------------- | :----------------------------------------------------------: | :----------------------------------------------------------- |
| **Quality Control (PCA)**  |       <img src="figures/001_QC_PCA.png" width="300"/>        | High sample reproducibility; clear separation between NC and siRNA groups. |
| **DEG Identification**     |       <img src="figures/002_Volcano.png" width="300"/>       | Volcano plot displaying the distribution of significantly altered genes. |
| **GO-BP (Down-regulated)** |    <img src="figures/004-1_GO_BP_down.png" width="300"/>     | Down-regulated genes are enriched in immune and cytokine-related processes. |
| **KEGG Pathway Mapping**   | <img src="figures/008-1_hsa04060.pathview.png" width="350"/> | Down-regulated genes (green boxes) highlighted on the hsa04060 pathway map. |

---



## Directory Structure

```text
.
├── README.md                			 # Project documentation
├── EXOSC10-AS1.Rproj         		     # RStudio project file
├── data/                    			 # Raw count matrices and exported DEG tables
│   ├── GSE289213_Genes_expression.xlsx
│   ├── DEG_results_siRNA1_vs_NC.csv
│   └── GO_BP_down_results.csv
├── figures/                 			  # High-resolution exported figures
│   ├── 001_QC_PCA.png
│   ├── 002_Volcano.png
│   ├── 003-1_up_Top10_heatmap.png
│   ├── ...
│   └── 008-1_hsa04060.pathview.png
└── scripts/                			  # Modularized R scripts
    ├── 01_prepare_and_QC.R  			  # Data preprocessing and PCA quality control
    ├── 02_volcano_and_top10.R		  	  # DESeq2 calculation and heatmap generation
    ├── 03_GO_enrichment.R   			  # Sub-dimensional GO enrichment analysis
    └── 04_KEGG_and_pathview.R			  # KEGG enrichment and Pathview mapping
```



## Requirements & Reproducibility (Quick Start)

### 1. Dependency Installation

Automated package management via `pak` is integrated within the scripts. Alternatively, dependencies can be manually installed:



```
# R

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "pathview", "pheatmap"))
install.packages(c("tidyverse", "readxl"))
```

### 2. Execution Pipeline

Create an R Project in the root directory and execute the 4 scripts in sequence to reproduce all outputs:

1. `01_prepare_and_QC.R`
2. `02_volcano_and_top10.R`
3. `03_GO_enrichment.R`
4. `04_KEGG_and_pathview.R`

## Technical Challenges & Troubleshooting

During the pipeline implementation, key data formatting, visualization, and API connection challenges were identified and addressed:

### 1. Integer Data Validation and DESeq2 Fault Tolerance

- **Issue**: `DESeq2` strictly requires raw count matrices to consist of non-negative integers. Floating-point values or incomplete type conversions trigger fatal errors in `DESeqDataSetFromMatrix()`.

- **Resolution**: Implemented automated validation and type casting prior to constructing the `dds` object:

  

  ```
  # R
  stopifnot(all(raw_df[, -1] %% 1 == 0))
  mode(counts_matrix) <- "integer"
  ```

### 2. Decoupling Volcano Plot Text Overlaps for Better Visualization

- **Issue**: Direct text labeling (`ggrepel`) of top DEGs on a single volcano plot led to severe label overlaps, cluttering data points and diminishing readability.
- **Resolution**: Applied a visualization decoupling strategy. Kept an uncluttered standard volcano plot for global distribution, while isolating top 10 up-/down-regulated genes into separate clustered heatmaps for in-depth analysis.

### 3. Functional Masking in Combined Up/Down DEG Analysis

- **Issue**: Initial GO-BP analysis combining both up- and down-regulated DEGs resulted in functional interference between upregulated (cellular compensatory) and downregulated (immune silencing) signals, obscuring the primary mechanism.
- **Resolution**: Adopted a **split enrichment strategy**, segregating genes into downregulated (193) and upregulated (90) subsets. This revealed that `EXOSC10-AS1` knockdown directly shuts down cytokine production and vesicle lumen structures.

### 4. Absence of Significant Enrichment in Upregulated Gene Sets

- **Issue**: GO-BP enrichment on upregulated genes (90 genes) yielded zero statistically significant terms under standard thresholds (`pAdjustMethod = "BH"`, `p.adjust < 0.05`), preventing dotplot generation.
- **Resolution**: Code and logic checks confirmed this was driven by data characteristics—upregulated genes were few and functionally dispersed, failing to form significant clusters. To maintain scientific rigor, strict statistical thresholds were preserved rather than artificially relaxed, focusing reporting on "Global GO-BP overview" and "Downregulated GO-BP/CC/MF deep dive".

### 5. Network Latency & API Timeout in KEGG Online Queries

- **Issue**: `clusterProfiler::enrichKEGG()` and `pathview` fetch real-time data from official KEGG API servers. Unstable network connectivity frequently caused timeout or SSL handshake errors.
- **Resolution**: Included network diagnostic checks. When encountering connection timeouts, configuring proxy settings in RStudio or switching VPN connections to TUN/Global mode ensures stable REST API access to KEGG servers.

### 6. Ensembl ID Mapping Deficits in KEGG Annotations

- **Issue**: `enrichKEGG()` cannot natively parse Ensembl IDs (`ENSG...`), and direct conversions risk data loss due to unannotated IDs in the NCBI database.
- **Resolution**: Utilized `clusterProfiler::bitr()` for precise ID mapping with `na.omit()` filtering, validating output conversion coverage to guarantee downstream KEGG enrichment integrity.


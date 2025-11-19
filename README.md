# DGE-Analysis-of-GSE166044-Using-DESeq2-with-PCA-Heatmap-and-Volcano-Plot-Visualization
DGE Analysis of GSE166044 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization
# 🧬 DESeq2 DGE Pipeline: Breast Tissue Prior to Cancer Diagnosis (GSE166044)

This R script automates a robust and complex bioinformatics workflow for **Differential Gene Expression (DGE)** analysis of the **GSE166044** dataset. This study examines gene expression in breast tissues **prior to cancer diagnosis**, comparing tissues deemed **Susceptible** to those classified as **Normal**.

The pipeline is specifically designed to handle the challenging data structure of this GEO entry (a TAR archive of individual count files) and provides comprehensive quality control (QC) and visualization using **`DESeq2`**.

## 🚀 Key Features

* **Complex Data Handling:** Automatically downloads, extracts (`untar`), and processes the **`GSE166044_RAW.tar`** archive containing individual count files.
* **Count Matrix Assembly:** Correctly iterates through all sample files, filters out non-gene summary lines, and assembles the final count matrix.
* **DESeq2 Analysis:** Implements the standard `DESeq2` workflow for robust DGE testing.
* **QC Transformation:** Applies **Variance Stabilizing Transformation (VST)** for accurate clustering and dimension reduction plots.
* **Integrated Visualization:** Generates three essential plots for interpreting the results: **PCA**, **Heatmap**, and **Volcano Plot**.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE166044 | Transcriptome analysis of breast tissues (Susceptible vs. Normal) prior to cancer diagnosis. |
| **Data Source** | GEO Supplementary Files (`.tar` archive) | Handles the complex file packaging unique to this dataset. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for RNA-Seq count data. |
| **Comparison** | Susceptible Tissue vs. Normal Tissue | Identifies early transcriptional changes in at-risk breast tissue. |
| **Significance** | $\text{pCutoff} = 0.05$, $\text{FCcutoff} = 1.5$ | Used for highlighting significant findings in the Volcano Plot. |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically checks for and installs the necessary Bioconductor and CRAN packages:
* `GEOquery` (For data download)
* `DESeq2` (For DGE analysis)
* `pheatmap` (For Heatmap visualization)
* `EnhancedVolcano` (For Volcano Plot visualization)
* `ggplot2` (For PCA visualization)
* `utils` (For `untar()` function)

### ⚙️ Execution

1.  **Download** the `DGE Analysis of GSE166044 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization.R` file.
2.  **Optional:** Modify the `output_dir` variable (Step 2b) to your preferred saving location.
    ```R
    output_dir <- "D:/DOWNLOADS" # Change this path
    ```
3.  **Execute** the script in your R environment:
    ```R
    source("DGE Analysis of GSE166044 Using DESeq2 with PCA, Heatmap, and Volcano Plot Visualization.R")
    ```
    *Note: The script will automatically download and extract the required supplementary files from GEO into a subdirectory named `GSE166044`.*

---

## 📁 Output Files (3 Plots)

The script automatically saves the following files to the specified `output_dir` (default: `D:/DOWNLOADS`).

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `GSE166044_PCA_plot.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating global clustering and separation of Susceptible vs. Normal samples. |
| `GSE166044_Heatmap_Top50.png` | QC | **Heatmap of the Top 50 Most Variable Genes** (based on VST-transformed data) to visualize sample grouping quality. |
| `GSE166044_Volcano_plot.png` | Results | **Volcano Plot** generated using `EnhancedVolcano`, showing the $\log_2 \text{Fold Change}$ vs. $P_{\text{value}}$, highlighting significant and highly changed genes. |

# Single-Cell RNA-Seq & Cell-Cell Communication Pipeline (NSCLC)

This repository contains a comprehensive, end-to-end bioinformatics workflow written in **R** for analyzing single-cell RNA sequencing (scRNA-seq) data. The project processes a high-throughput **Non-Small Cell Lung Cancer (NSCLC)** dataset using **Seurat** for quality control, scaling, and dimensionality reduction, and subsequently transitions into **CellChat** to model and visualize complex intercellular communication networks.

---

## Repository Workflow Structure

### Part 1: Seurat Preprocessing & Dimensionality Reduction
* **Data Ingestion:** Loads raw 10x Genomics multi-modal `.h5` expression matrices (`Read10X_h5`) and isolates the primary **Gene Expression** modality matrix.
* **Quality Control (QC):** Filters cells based on feature constraints (`nFeature_RNA > 200 & < 2500`) and mitochondrial gene content (`percent.mt < 5%`).
* **Normalization & Variable Features:** Applies standard log-normalization and extracts the top 2,000 highly variable features using the `vst` selection method.
* **Scaling & PCA:** Centers and scales data across all features, followed by linear dimensionality reduction via Principal Component Analysis (PCA).
* **Clustering & UMAP:** Constructs a K-nearest neighbor graph (using the top 15 principal components), identifies clusters across multiple resolutions (highlighting resolution `0.1`), and projects cells into a 2D UMAP space.

### Part 2: CellChat Communication Inference
* **Object Initialization:** Extracts normalized expression data and cluster identities from the Seurat object to build a structured CellChat framework.
* **Database Integration & Filtering:** Utilizes the human ligand-receptor interaction database (`CellChatDB.human`), focusing specifically on **Secreted Signaling** pathways.
* **Probability & Pathway Computation:** Computes communication probabilities, filters low-cell-count interactions (`min.cells = 10`), and aggregates networks at both individual ligand-receptor pair and system signaling pathway levels.

### Part 3: Downstream Network Visualization
The pipeline generates multi-faceted visualizations to interpret intercellular signals, including:
* **Circular Network Plots:** Showcasing overall interaction counts and weights/strengths across cell groups, as well as sender-specific pathway distributions.
* **Contribution Analyses:** Highlighting the fractional contribution of individual ligand-receptor pairs within specific signaling pathways (e.g., `MIF`, `SPP1`, `TNF`).
* **Chord Diagrams & Heatmaps:** Providing clear mappings of directional signaling dynamics and pathway heatmaps (e.g., `Reds` color palette for `TNF`).
* **Violin & Gene Expression Plots:** Verifying individual gene expression levels driving specific ligand-receptor interactions.

---

## Prerequisites & Dependencies

To execute this script successfully, ensure your R environment has the following packages installed:

```R
# Install core packages if not already present
install.packages("tidyverse")
install.packages("Seurat")

# Install CellChat from GitHub (if required)
# devtools::install_github("sqjin/CellChat")

library(Seurat)
library(tidyverse)
library(CellChat)

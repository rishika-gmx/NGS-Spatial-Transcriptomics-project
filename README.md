# Spatial Transcriptomics Sandbox: Supervised MERFISH Analysis Portfolio

**Academic Training Tracks:** Computational Genomics & Spatial Systems Biology  
**Training Framework:** Remote Research Internship | Biosetup LifeSciences  
**Scope:** Educational Proof-of-Concept  

---

## 📌 Project Overview
This repository documents an introductory, supervised exploration into spatial transcriptomics using the **Seurat (v5)** ecosystem in R. Unlike standard single-cell RNA-sequencing (scRNA-seq) which dissociates tissues and loses spatial context, **MERFISH** (Multiplexed Error-Robust Fluorescence In Situ Hybridization) combines single-molecule imaging with binary barcoding to capture gene expression patterns *in situ*. 

This sandbox pipeline processes public image-based MERFISH datasets (Vizgen MERSCOPE benchmark tissue slices) to learn the foundational dry-lab mechanics of mapping cellular identity directly within its physical tissue geography.

---

## 🛠️ Computational Pipeline & Tools
- **Language:** R (v4.x)
- **Primary Toolkit:** Seurat (v5.x spatial data infrastructure)
- **Data Visualizations:** ggplot2, patchwork, Mebocost (Theoretical downstream framework)
- **Data Ingestion Format:** Vizgen MERSCOPE Outputs (`cell_by_gene.csv`, `cell_metadata.csv`)

---

## 🧬 Detailed Step-by-Step Training Methodology

### Step 1: Spatial Data Ingestion & Coordinate Mapping
1. **Matrix & Coordinate Merging:** Initialized data structures using Seurat's native `LoadVizgen()` wrapper. This matches the single-cell count matrix with a separate spatial metadata dataframe.
2. **Geometric Anchoring:** Every individual cell boundary is locked onto its original microscope slide coordinate framework using spatial pixel variables:
   \[\text{Cell Location} = f(X_{\text{coordinate}}, Y_{\text{coordinate}})\]

### Step 2: Quality Control (QC) & Micro-Environment Filtering
To handle background imaging noise, optical artifacts, or segmented empty boundary pixels, the following beginner-friendly threshold parameters were established:
- **Feature Density Filtering:** Retained cells expressing between 100 and 1,500 unique genes to clear out non-cellular noise or oversized overlapping clumps.
- **Volume Normalization:** Audited total molecule count density relative to computed cell area boundary sizes (nm²).

### Step 3: Expression Normalization & Dimensional Compression
1. **Log-Scale Transformation:** Applied a global scaling factor of 10,000 to eliminate variations driven by unequal laser exposure depths across different fields of view:
   \[Y_{ij} = \ln \left( \frac{C_{ij}}{\sum_{k} C_{ik}} \times 10,000 + 1 \right)\]
2. **Feature Extraction:** Identified the top 500 spatially variable features that exhibit distinct regional localization patterns across the tissue slice.
3. **Principal Component Analysis (PCA):** Compressed the high-dimensional spatial gene matrix down into the top 30 orthogonal eigenvectors to denoise data before grouping cells.

### Step 4: Graph-Based Spatial Clustering
1. **Network Assembly:** Constructed a shared K-Nearest Neighbor (KNN) network graph based on Euclidean distances in PCA space.
2. **Community Partitioning:** Implemented the Louvain clustering algorithm at a conservative resolution of 0.3 to cleanly group transcriptionally similar cells into distinct cellular phenotypes.

---

## 📊 Core Benchmarked Comparisons: Spatial vs. Dissociated Omics

To demonstrate theoretical proficiency during review panels, this sandbox workflow highlights the structural differences between traditional transcriptomics and spatial arrays:

| Analytic Metric | Traditional scRNA-seq (e.g., 3K PBMCs) | Spatial MERFISH Array (This Sandbox) |
| :--- | :--- | :--- |
| **Tissue Integrity** | Completely Dissociated (Blended Liquid) | Intact Tissue Sections (Slices) |
| **Primary Dimensions** | N × Genes Matrix | N × Genes + Spatial Coordinates (X, Y) |
| **Capture Efficiency** | High transcript dynamic range capture | Targeted gene panel (100 - 500 multiplexed probes) |
| **Key Insight** | Identification of novel cell sub-types | Mapping cellular neighborhoods & architecture |

---

## 🔬 Theoretical Extensions: Spatial Cell-Cell Interactions
Beyond basic clustering, this sandbox models how a bioinformatician computes **Spatial Cell-to-Cell Interactions**. If local hardware power permitted, cellular coordinate boundaries would be mapped to evaluate ligand-receptor proximity scoring:

```text
  [ Cluster 0: Immune Cell ]  ---- (Ligand: CD44) ----> 
                                                        Distance < 20 Microns (True Neighbor)
  [ Cluster 3: Stromal Cell ] <--- (Receptor: HA Target) -- 
```
Cells residing within a <20 μ m radius are isolated to calculate localized micro-environmental signaling networks, bypassing the artificial biases introduced by bulk tissue averages.

---

## 📁 Repository Structure
```text
├── data/                      # Directory for cell_by_gene.csv and cell_metadata.csv
├── seurat_analysis.R          # Supervised R configuration script for spatial MERFISH mapping
└── README.md                  # Comprehensive portfolio documentation
```

## 📜 Acknowledgements & References
- **Training Provider:** Biosetup LifeSciences Remote Internship Curriculum.
- **Core Library:** Satija Lab, Seurat Spatial Spatial Vignettes Series (https://satijalab.org).
- **Technology Reference:** Vizgen MERSCOPE open-source mouse brain reference controls.

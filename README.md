# NGS & Spatial Transcriptomics Analysis using Seurat

**Undergraduate Bioinformatics Internship Project**

Educational Proof-of-Concept | R • Seurat • Linux • MERFISH
---

## 📌 Project Overview
This repository documents an introductory bioinformatics project completed during a computational biology internship, focusing on the analysis of publicly available MERFISH spatial transcriptomics datasets using R and the Seurat package.

The project introduces a standard single-cell and spatial transcriptomics workflow, including quality control, normalization, dimensionality reduction, clustering, marker identification, and spatial visualization. It was completed as a supervised learning project to develop practical experience with transcriptomic data analysis and reproducible computational workflows.

---

## 🎯 Project Scope
This project demonstrates:

- Linux-based bioinformatics workflow
- R programming
- Seurat analysis
- Quality control
- PCA
- Clustering
- Marker gene identification
- Spatial transcriptomics visualization

This repository represents an educational implementation completed under supervision during a computational biology internship.

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


## 🔬 Analysis Workflow
```
Public MERFISH Dataset
            │
            ▼
      Data Import
            │
            ▼
     Quality Control
            │
            ▼
      Normalization
            │
            ▼
Dimensional Reduction (PCA)
            │
            ▼
    Nearest Neighbours
            │
            ▼
        Clustering
            │
            ▼
 Marker Gene Detection
            │
            ▼
 Spatial Visualization
```

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

## ✔ Implementation Status

```
Spatial Transcriptomics Workflow
│
├── Data import                          ✅ Completed
├── Quality control                      ✅ Completed
├── Normalization                        ✅ Completed
├── PCA                                  ✅ Completed
├── Nearest neighbor graph               ✅ Completed
├── Clustering                           ✅ Completed
├── Marker gene analysis                 ✅ Completed
├── Spatial visualization                ✅ Completed
└── Workflow documentation               ✅ Completed
```

## 📊 Expected Outputs

The analysis workflow generates:

- Quality control metrics
- PCA visualization
- Cluster assignments
- UMAP projection
- Marker gene tables
- Spatial feature plots
---

## 📁 Repository Structure

```
NGS-Spatial-Transcriptomics-project
├── README.md
├── scripts/
│   └── Seurat_analysis.R
├── docs/
├── figures/
├── results/
├── LICENSE
```

## 📜 Acknowledgements & References
- **Training Provider:** Biosetup LifeSciences Remote Internship Curriculum.
- **Core Library:** Satija Lab, Seurat Spatial Spatial Vignettes Series (https://satijalab.org).
- **Technology Reference:** Vizgen MERSCOPE open-source mouse brain reference controls.

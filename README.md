
````markdown
# NGS & Spatial Transcriptomics Analysis using Seurat

**Undergraduate Bioinformatics Internship Project**

Reproducible Portfolio Analysis | R • Seurat • Linux • MERFISH

---

## 📌 Project Overview

This repository documents a supervised introductory bioinformatics project completed during a computational biology internship. The project focuses on the analysis of a publicly available MERFISH spatial transcriptomics dataset from mouse brain tissue using **R and Seurat**.

The workflow was designed to provide practical experience with transcriptomic data processing, quality control, dimensionality reduction, graph-based clustering, marker-gene identification, spatial coordinate analysis, and biological interpretation of transcriptionally distinct cell populations.

The project was completed as a supervised learning exercise and is intended to demonstrate practical familiarity with a reproducible computational biology workflow rather than independent research.

---

## 🎯 Project Objectives

The analysis demonstrates practical experience with:

- R-based transcriptomic data analysis
- Linux-based computational workflows
- Seurat / Seurat v5
- MERFISH spatial transcriptomics data
- Quality control and cell filtering
- Expression normalization
- Highly variable gene identification
- Principal Component Analysis (PCA)
- Nearest-neighbour graph construction
- Graph-based clustering
- UMAP dimensionality reduction
- Marker-gene identification
- Marker-expression visualization
- Spatial coordinate analysis
- Cluster-level biological interpretation
- Manual and provisional cell-type annotation
- Reproducible generation of analysis outputs

---

## 🛠️ Computational Environment & Tools

- **Programming Language:** R
- **Primary Analysis Framework:** Seurat / Seurat v5
- **Visualization:** ggplot2, patchwork, Seurat plotting functions
- **Data Type:** MERFISH spatial transcriptomics
- **Dataset:** Publicly available mouse brain MERFISH dataset
- **Operating Environment:** Linux / WSL
- **Development Environment:** Visual Studio Code with R integration

---

## 🧬 Dataset & Data Ingestion

The analysis uses a publicly available MERFISH mouse brain dataset containing:

- **483 genes**
- **83,546 cells**
- Spatial coordinates for the individual cells

The supplied RData file contains the expression matrix and corresponding spatial coordinate information.

The initial data structures were inspected and validated before constructing the Seurat object.

Spatial coordinates were subsequently associated with Seurat cluster identities to enable visualization of the transcriptional clusters across the original tissue coordinate space.

---

# 🔬 Analysis Workflow

The project follows the workflow below:

```text
Public MERFISH Mouse Brain Dataset
              │
              ▼
       Data Ingestion
              │
              ▼
       Seurat Object
              │
              ▼
     Quality Control
              │
              ▼
       Normalization
              │
              ▼
 Highly Variable Features
              │
              ▼
             PCA
              │
              ▼
   Nearest-Neighbour Graph
              │
              ▼
        Clustering
              │
              ▼
            UMAP
              │
              ▼
    Marker Gene Detection
              │
              ▼
      Marker Heatmaps
              │
              ▼
   Spatial Cluster Mapping
              │
              ▼
 Cluster-Level Interpretation
              │
              ▼
 Manual / Provisional
 Cell-Type Annotation
````

---

## 1. Data Ingestion & Seurat Object Construction

The publicly available MERFISH dataset was loaded from an RData file containing:

* The gene-expression count matrix
* Cell-level spatial coordinates

The expression matrix was used to construct a Seurat object while preserving the correspondence between cells and their spatial coordinates.

The initial dataset contained:

**483 genes × 83,546 cells**

The spatial coordinate information was retained separately and later incorporated into the spatial cluster visualization workflow.

---

## 2. Quality Control

Quality-control metrics were generated from the expression matrix before downstream analysis.

The primary filtering criterion used in the workflow was based on the number of detected genes per cell (`nFeature_RNA`).

Cells with:

```text
nFeature_RNA < 20
```

were excluded from downstream analysis.

After filtering:

```text
Cells retained: 83,407
```

The QC stage was used to remove cells with extremely low detected-gene counts while retaining the majority of the dataset for subsequent analysis.

QC-related summary information was saved for reproducibility.

---

## 3. Expression Normalization

Following quality control, expression values were normalized using Seurat's:

```r
NormalizeData(
    normalization.method = "LogNormalize",
    scale.factor = 10000
)
```

Log-normalization converts raw expression counts into normalized expression values while accounting for differences in sequencing/read depth between cells.

---

## 4. Highly Variable Feature Identification

Highly variable genes were identified using Seurat's variance-stabilizing transformation (`vst`) approach.

The analysis retained:

```text
200 highly variable features
```

These features were subsequently used for dimensionality reduction and downstream transcriptional analysis.

---

## 5. Data Scaling

The expression data were scaled before PCA so that genes could be compared on a standardized scale during dimensionality reduction.

Scaling was performed using Seurat's `ScaleData()` workflow.

---

## 6. Principal Component Analysis (PCA)

Principal Component Analysis was performed on the selected highly variable features.

The analysis generated:

```text
30 principal components
```

PCA was used to reduce the dimensionality of the expression matrix while preserving major sources of transcriptional variation.

The resulting principal components formed the basis for downstream neighbourhood graph construction, clustering, and UMAP visualization.

---

## 7. Nearest-Neighbour Graph Construction

A nearest-neighbour graph was constructed from the PCA representation of the cells.

The graph connects transcriptionally similar cells based on their positions in reduced-dimensional expression space.

This graph was subsequently used for Seurat's graph-based clustering algorithm.

---

## 8. Graph-Based Clustering

Seurat's graph-based clustering workflow was applied to identify transcriptionally similar groups of cells.

The clustering analysis used:

```text
Resolution: 0.5
```

The final analysis produced:

```text
20 transcriptional clusters
```

The resulting cluster assignments were used throughout the remainder of the analysis.

Cluster sizes were summarized and saved to:

```text
results/cluster_summary.csv
```

---

## 9. UMAP Dimensionality Reduction

UMAP was performed using the selected principal components to generate a two-dimensional representation of the transcriptional relationships between cells.

The current workflow uses:

```text
PCs 1–15
```

The resulting visualization provides an overview of the 20 transcriptional clusters and their relative organization in reduced-dimensional expression space.

Output:

```text
figures/umap_clusters.png
```

---

## 10. Marker Gene Identification

Cluster-specific marker genes were identified using Seurat's marker-detection workflow.

Positive marker genes were examined to determine which genes were enriched within individual transcriptional clusters.

The marker analysis produced cluster-level marker tables, including:

```text
results/marker_genes_all_clusters.csv
results/top10_markers_per_cluster.csv
results/top3_markers_per_cluster.csv
results/cluster_markers.csv
```

These tables were subsequently used for visualization and biological interpretation.

---

## 11. Marker Gene Heatmaps

Marker-gene expression was visualized across the identified clusters using heatmaps.

Two complementary heatmaps were generated.

### Marker Expression Heatmap

The cluster-level marker heatmap displays expression patterns of selected marker genes across individual cells grouped by cluster.

Output:

```text
figures/marker_heatmap.png
```

### Average Marker Expression Heatmap

A second heatmap summarizes average marker-gene expression across each cluster.

This provides a cleaner cluster-level view of marker enrichment and is useful for comparing transcriptional signatures between clusters.

Output:

```text
figures/marker_heatmap_average_expression.png
```

The heatmaps were used as complementary visualization tools rather than as independent clustering methods.

---

## 12. Spatial Coordinate Analysis

The MERFISH dataset contains spatial coordinates corresponding to individual cells.

Rather than relying on a native tissue-image object, the project uses the supplied spatial coordinate information to map Seurat cluster identities back onto the original tissue coordinate system.

The spatial coordinates were extracted and saved as:

```text
results/spatial_coordinates.csv
```

Cluster identities were then combined with the spatial coordinates to visualize where transcriptionally distinct populations occur within the tissue.

Output:

```text
figures/spatial_clusters.png
```

This provides a spatial interpretation of the transcriptional clusters identified through the Seurat workflow.

---

# 🧠 Cluster-Level Biological Interpretation

Following marker-gene analysis, clusters were manually reviewed using their marker profiles.

The current annotation workflow identifies broad cell-type identities where the marker evidence is sufficiently strong while retaining **provisional or unresolved labels** for clusters where the available marker profile does not support a confident subtype assignment.

The resulting annotation files include:

```text
results/preliminary_cluster_annotations.csv
results/cluster_annotations_manual.csv
```

The final annotations include broad populations such as:

* Excitatory / glutamatergic neurons
* GABAergic / inhibitory neurons
* Oligodendrocyte-lineage cells
* Oligodendrocyte precursor cells (OPCs)
* Astrocytes
* Endothelial cells
* Microglia
* Vascular mural cells / pericytes
* Vascular smooth muscle cells
* Cholinergic neurons
* Other specialized or unresolved neuronal populations

These annotations represent **manual biological interpretation of the computational clusters** and should not be interpreted as experimentally validated cell identities.

---

# 📊 Current Cluster Summary

The final clustering workflow identified **20 transcriptional clusters**.

| Cluster | Cell Count |
| ------: | ---------: |
|       0 |     12,258 |
|       1 |      8,821 |
|       2 |      8,545 |
|       3 |      8,219 |
|       4 |      7,477 |
|       5 |      6,774 |
|       6 |      5,393 |
|       7 |      4,452 |
|       8 |      4,139 |
|       9 |      3,457 |
|      10 |      3,328 |
|      11 |      2,426 |
|      12 |      2,236 |
|      13 |      1,993 |
|      14 |      1,057 |
|      15 |        765 |
|      16 |        735 |
|      17 |        679 |
|      18 |        475 |
|      19 |        178 |

The complete machine-readable summary is available in:

```text
results/cluster_summary.csv
```

---

# 📁 Repository Structure

```text
NGS-Spatial-Transcriptomics-project/
│
├── data/
│   └── merfish_mouse_brain/
│       └── Vizgen_Merfish_count_location.RData
│
├── figures/
│   ├── marker_heatmap.png
│   ├── marker_heatmap_average_expression.png
│   ├── spatial_clusters.png
│   └── umap_clusters.png
│
├── results/
│   ├── cluster_annotations_manual.csv
│   ├── cluster_markers.csv
│   ├── cluster_summary.csv
│   ├── curated_heatmap_markers.csv
│   ├── marker_genes_all_clusters.csv
│   ├── preliminary_cluster_annotations.csv
│   ├── qc_summary.csv
│   ├── sessionInfo.txt
│   ├── spatial_coordinates.csv
│   ├── top10_markers_per_cluster.csv
│   ├── top3_markers_per_cluster.csv
│   ├── merfish_seurat_processed.rds
│   └── merfish_seurat_annotated.rds
│
├── .gitignore
├── README.md
├── Rplots.pdf
├── Seurat_analysis.R
└── LICENSE
```

---

# 📦 Generated Outputs

The analysis generates the following major outputs.

### Figures

```text
figures/
├── marker_heatmap.png
├── marker_heatmap_average_expression.png
├── spatial_clusters.png
└── umap_clusters.png
```

### Cluster and marker results

```text
results/
├── cluster_summary.csv
├── cluster_markers.csv
├── marker_genes_all_clusters.csv
├── top3_markers_per_cluster.csv
├── top10_markers_per_cluster.csv
└── curated_heatmap_markers.csv
```

### Annotation results

```text
results/
├── preliminary_cluster_annotations.csv
└── cluster_annotations_manual.csv
```

### Quality control and spatial information

```text
results/
├── qc_summary.csv
└── spatial_coordinates.csv
```

### Reproducibility and processed objects

```text
results/
├── sessionInfo.txt
├── merfish_seurat_processed.rds
└── merfish_seurat_annotated.rds
```

---

# ✔ Implementation Status

```text
MERFISH Spatial Transcriptomics Workflow
│
├── Public dataset ingestion                 ✅ Completed
├── Seurat object construction               ✅ Completed
├── Quality control                          ✅ Completed
├── Expression normalization                 ✅ Completed
├── Highly variable feature selection        ✅ Completed
├── PCA                                      ✅ Completed
├── Nearest-neighbour graph                  ✅ Completed
├── Graph-based clustering                   ✅ Completed
├── UMAP visualization                       ✅ Completed
├── Marker-gene identification               ✅ Completed
├── Marker heatmap                           ✅ Completed
├── Average-expression heatmap               ✅ Completed
├── Spatial coordinate mapping               ✅ Completed
├── Spatial cluster visualization            ✅ Completed
├── Cluster-level marker interpretation      ✅ Completed
├── Manual/provisional annotation            ✅ Completed
└── Reproducibility documentation            ✅ Completed
```

---

# 🔬 Future Extensions

The current project focuses on the core transcriptomic and spatial-clustering workflow.

Potential extensions include:

* More detailed neuronal subtype annotation
* Additional marker validation using external cell-type references
* Differential expression analysis between selected populations
* Spatial neighbourhood analysis
* Quantitative cell-cell proximity analysis
* Ligand-receptor interaction analysis
* Spatially variable gene analysis
* Integration with additional transcriptomic datasets

These analyses are **not part of the current completed workflow** and are listed only as potential future directions.

---

# 📚 Reproducibility

The complete analysis is documented in:

```text
Seurat_analysis.R
```

The processed Seurat objects and intermediate result tables are provided in the `results/` directory.

R session and package information are recorded in:

```text
results/sessionInfo.txt
```

The analysis was performed in a Linux/WSL environment using R and Seurat.

---

# 📜 Acknowledgements & References

* **Training Provider:** Biosetup LifeSciences Remote Internship Curriculum.
* **Core Analysis Framework:** Seurat / Satija Lab.
* **Spatial Transcriptomics Technology:** Vizgen MERFISH / MERSCOPE.
* **Dataset:** Publicly available MERFISH mouse brain dataset.

This project was completed as a supervised computational biology learning exercise and is presented for educational, portfolio, and reproducibility purposes.


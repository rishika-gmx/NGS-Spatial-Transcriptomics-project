# Spatial Transcriptomics Analysis using MERFISH Data
#
# This script documents a supervised spatial transcriptomics workflow
# developed during a computational biology internship.
#
# The workflow is reproduced using a publicly available MERFISH dataset
# for portfolio and reproducibility purposes.
#
# Target Toolsets:
#   - R
#   - Seurat v5
#   - ggplot2
#   - patchwork
#   - Vizgen/MERSCOPE spatial data
#
# -------------------------------------------------------------------------
# 1. Load Required Libraries
# -------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(patchwork)

cat("\n")
cat("=========================================================================\n")
cat("--- MERFISH Spatial Transcriptomics Analysis Pipeline -------------------\n")
cat("=========================================================================\n\n")

# -------------------------------------------------------------------------
# 2. Define Input and Output Paths
# -------------------------------------------------------------------------

# Directory containing the publicly available MERFISH dataset.
#
# Expected location:
# data/merfish_mouse_brain/
#
# The directory should contain the files required by LoadVizgen().

data_directory <- "data/merfish_mouse_brain/"

# Create directories for reproducible outputs.
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# 3. Load MERFISH Spatial Transcriptomics Data
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 1: Spatial Data Ingestion --------------------------------------\n")
cat("=========================================================================\n")

if (!dir.exists(data_directory)) {

  stop(
    paste0(
      "\nMERFISH dataset not found.\n\n",
      "Expected dataset directory:\n",
      data_directory,
      "\n\n",
      "Please place the publicly available MERFISH dataset in this location ",
      "before running the analysis."
    )
  )

}

cat("Loading MERFISH spatial transcriptomics data...\n")

merfish_data <- LoadVizgen(
  data.dir = data_directory
)

cat(
  "Successfully loaded MERFISH dataset with ",
  ncol(merfish_data),
  " cells and ",
  nrow(merfish_data),
  " features.\n\n"
)

# -------------------------------------------------------------------------
# 4. Quality Control
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 2: Quality Control ---------------------------------------------\n")
cat("=========================================================================\n")

cat("Generating baseline quality-control metrics...\n")

# Visualize the number of detected features per cell.
qc_plot <- VlnPlot(
  merfish_data,
  features = "nFeature_RNA",
  pt.size = 0
) +
  ggtitle("MERFISH Quality Control: Detected Features per Cell") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = "figures/qc_nFeature_RNA.png",
  plot = qc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Apply the filtering thresholds used in the original workflow.
# Cells with fewer than 100 detected features or more than 1500 detected
# features are excluded from downstream analysis.

merfish_data <- subset(
  merfish_data,
  subset = nFeature_RNA > 100 & nFeature_RNA < 1500
)

cat(
  "QC filtering complete.\n",
  "Cells retained: ",
  ncol(merfish_data),
  "\n\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 5. Normalization
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 3: Expression Normalization ------------------------------------\n")
cat("=========================================================================\n")

cat("Applying LogNormalize...\n")

merfish_data <- NormalizeData(
  merfish_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# -------------------------------------------------------------------------
# 6. Identify Highly Variable Features
# -------------------------------------------------------------------------

cat("Identifying highly variable genes...\n")

merfish_data <- FindVariableFeatures(
  merfish_data,
  selection.method = "vst",
  nfeatures = 500
)

cat(
  "Number of highly variable features identified: ",
  length(VariableFeatures(merfish_data)),
  "\n\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 7. Scale Data
# -------------------------------------------------------------------------

cat("Scaling expression data...\n")

all_genes <- rownames(merfish_data)

merfish_data <- ScaleData(
  merfish_data,
  features = all_genes
)

# -------------------------------------------------------------------------
# 8. Principal Component Analysis
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 4: Principal Component Analysis -------------------------------\n")
cat("=========================================================================\n")

cat("Running PCA...\n")

merfish_data <- RunPCA(
  merfish_data,
  features = VariableFeatures(merfish_data),
  npcs = 30,
  verbose = FALSE
)

pca_plot <- DimPlot(
  merfish_data,
  reduction = "pca",
  label = TRUE
) +
  ggtitle("PCA of MERFISH Expression Profiles") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = "figures/pca_clusters.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 9. Construct Nearest-Neighbour Graph
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 5: Graph-Based Clustering --------------------------------------\n")
cat("=========================================================================\n")

cat("Constructing nearest-neighbour graph...\n")

merfish_data <- FindNeighbors(
  merfish_data,
  dims = 1:10
)

# -------------------------------------------------------------------------
# 10. Identify Cell Clusters
# -------------------------------------------------------------------------

cat("Identifying transcriptionally similar cell populations...\n")

merfish_data <- FindClusters(
  merfish_data,
  resolution = 0.3,
  verbose = FALSE
)

cat(
  "Number of identified clusters: ",
  length(unique(Idents(merfish_data))),
  "\n\n",
  sep = ""
)

# -------------------------------------------------------------------------
# 11. UMAP Dimensional Reduction
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 6: UMAP Visualization ------------------------------------------\n")
cat("=========================================================================\n")

cat("Running UMAP...\n")

merfish_data <- RunUMAP(
  merfish_data,
  dims = 1:10
)

umap_plot <- DimPlot(
  merfish_data,
  reduction = "umap",
  label = TRUE
) +
  ggtitle("UMAP of MERFISH Cell Clusters") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  filename = "figures/umap_clusters.png",
  plot = umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------------------------------------------------------
# 12. Spatial Visualization
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 7: Spatial Cluster Visualization -------------------------------\n")
cat("=========================================================================\n")

# LoadVizgen() stores spatial information when the supplied dataset contains
# the required spatial coordinate information.

if ("Vizgen" %in% names(merfish_data@images)) {

  cat("Generating spatial cluster map...\n")

  spatial_plot <- SpatialDimPlot(
    merfish_data,
    stroke = NA
  ) +
    ggtitle("Spatial Distribution of MERFISH Cell Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(
    filename = "figures/spatial_clusters.png",
    plot = spatial_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

} else {

  cat(
    "Spatial image information was not detected in the loaded object.\n",
    "Spatial visualization was therefore skipped.\n"
  )

}

# -------------------------------------------------------------------------
# 13. Identify Marker Genes
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 8: Marker Gene Identification ----------------------------------\n")
cat("=========================================================================\n")

cat("Identifying cluster marker genes...\n")

markers <- FindAllMarkers(
  merfish_data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save the complete marker-gene table.
write.csv(
  markers,
  file = "results/marker_genes.csv",
  row.names = FALSE
)

cat(
  "Marker-gene analysis complete.\n",
  "Results saved to results/marker_genes.csv\n\n"
)

# -------------------------------------------------------------------------
# 14. Cluster Summary
# -------------------------------------------------------------------------

cat("Generating cluster summary...\n")

cluster_counts <- as.data.frame(
  table(Idents(merfish_data))
)

colnames(cluster_counts) <- c(
  "Cluster",
  "Cell_Count"
)

write.csv(
  cluster_counts,
  file = "results/cluster_summary.csv",
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 15. Spatial Coordinate Audit
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Step 9: Spatial Coordinate Audit ------------------------------------\n")
cat("=========================================================================\n")

if ("Vizgen" %in% names(merfish_data@images)) {

  cat("Extracting spatial coordinate information...\n")

  spatial_coordinates <- GetTissueCoordinates(
    merfish_data,
    image = "Vizgen"
  )

  write.csv(
    spatial_coordinates,
    file = "results/spatial_coordinates.csv",
    row.names = FALSE
  )

  cat(
    "Spatial coordinates saved to results/spatial_coordinates.csv\n\n"
  )

} else {

  cat(
    "Spatial coordinate information was not available in the loaded object.\n\n"
  )

}

# -------------------------------------------------------------------------
# 16. Save Processed Seurat Object
# -------------------------------------------------------------------------

cat("Saving processed Seurat object...\n")

saveRDS(
  merfish_data,
  file = "results/merfish_processed.rds"
)


# -------------------------------------------------------------------------
# 17. Reproducibility Information
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Reproducibility Information -----------------------------------------\n")
cat("=========================================================================\n")

# Save R and package information used for this analysis.
capture.output(
  sessionInfo(),
  file = "results/sessionInfo.txt"
)

cat(
  "Session information saved to results/sessionInfo.txt\n\n"
)

# -------------------------------------------------------------------------
# 18. Analysis Complete
# -------------------------------------------------------------------------

cat("=========================================================================\n")
cat("--- Analysis Complete ---------------------------------------------------\n")
cat("=========================================================================\n")

cat("\nGenerated outputs:\n")
cat("  figures/qc_nFeature_RNA.png\n")
cat("  figures/pca_clusters.png\n")
cat("  figures/umap_clusters.png\n")
cat("  figures/spatial_clusters.png (if spatial information is available)\n")
cat("  results/marker_genes.csv\n")
cat("  results/cluster_summary.csv\n")
cat("  results/spatial_coordinates.csv (if available)\n")
cat("  results/merfish_processed.rds\n")
cat("  results/sessionInfo.txt\n\n")

cat("MERFISH spatial transcriptomics analysis completed successfully.\n")
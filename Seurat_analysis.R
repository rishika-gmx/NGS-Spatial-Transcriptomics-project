# ============================================================================
# Spatial Transcriptomics Analysis using MERFISH Data
# ============================================================================
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
#   - dplyr
#   - patchwork
#   - Vizgen/MERSCOPE spatial data
#
# Workflow:
#   1. Load libraries
#   2. Define paths
#   3. Load MERFISH data
#   4. Create Seurat object
#   5. Quality control
#   6. Normalization
#   7. Highly variable feature identification
#   8. Scaling
#   9. PCA
#  10. Nearest-neighbour graph
#  11. Graph-based clustering
#  12. UMAP
#  13. Spatial visualization
#  14. Marker-gene identification
#  15. Marker ranking and summary
#  16. Marker heatmaps
#  17. Manual cluster annotation
#  18. Spatial coordinate audit
#  19. Save annotated Seurat object
#  20. Reproducibility information
#
# ============================================================================


# ============================================================================
# 1. Load Required Libraries
# ============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

cat("\n")
cat("=========================================================================\n")
cat("--- MERFISH Spatial Transcriptomics Analysis Pipeline -------------------\n")
cat("=========================================================================\n\n")


# ============================================================================
# 2. Define Input and Output Paths
# ============================================================================

data_file <- "data/merfish_mouse_brain/Vizgen_Merfish_count_location.RData"

dir.create(
  "figures",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================================
# 3. Load MERFISH Spatial Transcriptomics Data
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 1: Spatial Data Ingestion --------------------------------------\n")
cat("=========================================================================\n")

if (!file.exists(data_file)) {

  stop(
    paste0(
      "\nMERFISH dataset not found.\n\n",
      "Expected file:\n",
      data_file,
      "\n\n",
      "Please check that the downloaded RData file is located in:\n",
      "data/merfish_mouse_brain/"
    )
  )

}

cat("Loading publicly available MERFISH dataset...\n")

load(data_file)


# Validate that the expected objects were loaded
if (!exists("raw_matrix") || !exists("location")) {

  stop(
    paste0(
      "\nThe RData file does not contain the expected objects ",
      "'raw_matrix' and 'location'."
    )
  )

}

cat(
  "Count matrix dimensions: ",
  nrow(raw_matrix),
  " genes x ",
  ncol(raw_matrix),
  " cells\n",
  sep = ""
)

cat(
  "Spatial coordinate dimensions: ",
  nrow(location),
  " cells x ",
  ncol(location),
  " coordinates\n\n",
  sep = ""
)


# ============================================================================
# 4. Create Seurat Object
# ============================================================================

cat("Creating Seurat object from MERFISH count matrix...\n")

merfish_data <- CreateSeuratObject(
  counts = raw_matrix,
  project = "MERFISH_Mouse_Brain",
  min.cells = 0,
  min.features = 0
)

cat(
  "Seurat object successfully created with ",
  ncol(merfish_data),
  " cells and ",
  nrow(merfish_data),
  " genes.\n\n",
  sep = ""
)


# ============================================================================
# 5. Quality Control
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 2: Quality Control ---------------------------------------------\n")
cat("=========================================================================\n")

cat("Generating baseline quality-control metrics...\n")


# Visualize the number of detected features per cell
qc_plot <- VlnPlot(
  merfish_data,
  features = "nFeature_RNA",
  pt.size = 0
) +
  ggtitle("MERFISH Quality Control: Detected Features per Cell") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "figures/qc_nFeature_RNA.png",
  plot = qc_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# Record QC statistics before filtering
qc_before <- data.frame(
  metric = "cells_before_filtering",
  value = ncol(merfish_data)
)


# Apply the filtering thresholds used in the workflow
#
# Cells with fewer than 100 detected features or more than 1500 detected
# features are excluded from downstream analysis.

merfish_data <- subset(
  merfish_data,
  subset = nFeature_RNA > 100 & nFeature_RNA < 1500
)


# Record QC statistics after filtering
qc_after <- data.frame(
  metric = "cells_after_filtering",
  value = ncol(merfish_data)
)

qc_summary <- rbind(
  qc_before,
  qc_after
)

write.csv(
  qc_summary,
  file = "results/qc_summary.csv",
  row.names = FALSE
)

cat(
  "QC filtering complete.\n",
  "Cells retained: ",
  ncol(merfish_data),
  "\n\n",
  sep = ""
)


# ============================================================================
# 6. Normalization
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 3: Expression Normalization ------------------------------------\n")
cat("=========================================================================\n")

cat("Applying LogNormalize...\n")

merfish_data <- NormalizeData(
  merfish_data,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ============================================================================
# 7. Identify Highly Variable Features
# ============================================================================

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


# ============================================================================
# 8. Scale Data
# ============================================================================

cat("Scaling expression data...\n")

all_genes <- rownames(merfish_data)

merfish_data <- ScaleData(
  merfish_data,
  features = all_genes
)


# ============================================================================
# 9. Principal Component Analysis
# ============================================================================

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
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "figures/pca_clusters.png",
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================================
# 10. Construct Nearest-Neighbour Graph
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 5: Graph-Based Clustering --------------------------------------\n")
cat("=========================================================================\n")

cat("Constructing nearest-neighbour graph...\n")

merfish_data <- FindNeighbors(
  merfish_data,
  dims = 1:10
)


# ============================================================================
# 11. Identify Cell Clusters
# ============================================================================

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


# ============================================================================
# 12. UMAP Dimensional Reduction
# ============================================================================

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
  ggtitle("MERFISH UMAP - Seurat Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "figures/umap_clusters.png",
  plot = umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)


# ============================================================================
# 13. Spatial Visualization
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 7: Spatial Cluster Visualization -------------------------------\n")
cat("=========================================================================\n")


# If Vizgen spatial information is available inside the Seurat object,
# generate the Seurat spatial plot.

if ("Vizgen" %in% names(merfish_data@images)) {

  cat("Generating spatial cluster map...\n")

  spatial_plot <- SpatialDimPlot(
    merfish_data,
    stroke = NA
  ) +
    ggtitle("MERFISH Spatial Distribution of Seurat Clusters") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(
    filename = "figures/spatial_clusters.png",
    plot = spatial_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

} else {

  cat(
    "Spatial image information was not detected in the Seurat object.\n",
    "Spatial visualization was therefore skipped.\n"
  )

}


# ============================================================================
# 14. Identify Marker Genes
# ============================================================================

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


# Save complete marker table
write.csv(
  markers,
  file = "results/marker_genes_all_clusters.csv",
  row.names = FALSE
)

cat(
  "Marker-gene analysis complete.\n",
  "Results saved to results/marker_genes_all_clusters.csv\n\n"
)


# ============================================================================
# 15. Marker Ranking and Cluster Summary
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 9: Marker Ranking and Cluster Summaries ------------------------\n")
cat("=========================================================================\n")


# --------------------------------------------------------------------------
# 15.1 Save simplified cluster marker table
# --------------------------------------------------------------------------

cluster_markers <- markers %>%
  select(
    cluster,
    gene,
    avg_log2FC,
    pct.1,
    pct.2,
    p_val_adj
  )

write.csv(
  cluster_markers,
  file = "results/cluster_markers.csv",
  row.names = FALSE
)


# --------------------------------------------------------------------------
# 15.2 Rank markers within each cluster
# --------------------------------------------------------------------------

markers_ranked <- markers %>%
  group_by(cluster) %>%
  arrange(
    desc(avg_log2FC),
    p_val_adj,
    .by_group = TRUE
  )


# --------------------------------------------------------------------------
# 15.3 Top 3 markers per cluster
# --------------------------------------------------------------------------

top3_markers <- markers_ranked %>%
  slice_head(n = 3) %>%
  ungroup()

write.csv(
  top3_markers,
  file = "results/top3_markers_per_cluster.csv",
  row.names = FALSE
)


# --------------------------------------------------------------------------
# 15.4 Top 10 markers per cluster
# --------------------------------------------------------------------------

top10_markers <- markers_ranked %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(
  top10_markers,
  file = "results/top10_markers_per_cluster.csv",
  row.names = FALSE
)


# --------------------------------------------------------------------------
# 15.5 Cluster summary
# --------------------------------------------------------------------------

cluster_counts <- as.data.frame(
  table(Idents(merfish_data))
)

colnames(cluster_counts) <- c(
  "cluster",
  "cell_count"
)

write.csv(
  cluster_counts,
  file = "results/cluster_summary.csv",
  row.names = FALSE
)

cat(
  "Marker rankings and cluster summary saved.\n\n"
)


# ============================================================================
# 16. Marker Gene Heatmaps
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 10: Marker Gene Heatmap Analysis -------------------------------\n")
cat("=========================================================================\n")


# --------------------------------------------------------------------------
# 16.1 Select curated marker genes
# --------------------------------------------------------------------------
#
# These genes were selected to represent major neuronal, glial and vascular
# populations identified during marker inspection.
#
# The same curated marker set is used for the heatmap outputs so that the
# figures remain directly comparable.

curated_markers <- c(
  "Adora2a",
  "Cckar",
  "Cx3cr1",
  "Drd1",
  "Drd2",
  "P2ry12",
  "S1pr5",
  "Gpr17",
  "Aldh1l1",
  "Olig1",
  "Gad1",
  "Slc17a6",
  "Slc17a7",
  "Aqp4",
  "Cldn5",
  "Myh11",
  "Kcnj8",
  "Abcc9",
  "Chat",
  "Slc32a1",
  "Gfap",
  "Csf1r",
  "Erbb3",
  "Erbb4",
  "Flt1",
  "Kdr",
  "Pdgfra",
  "Pdgfrb",
  "Gpr151"
)

curated_markers <- curated_markers[
  curated_markers %in% rownames(merfish_data)
]

write.csv(
  data.frame(gene = curated_markers),
  file = "results/curated_heatmap_markers.csv",
  row.names = FALSE
)


# --------------------------------------------------------------------------
# 16.2 Cell-level marker heatmap
# --------------------------------------------------------------------------

cat("Generating cell-level marker heatmap...\n")

marker_heatmap <- DoHeatmap(
  merfish_data,
  features = curated_markers,
  group.by = "seurat_clusters"
) +
  ggtitle("Marker Gene Expression Across Seurat Clusters") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "figures/marker_heatmap.png",
  plot = marker_heatmap,
  width = 14,
  height = 10,
  dpi = 300
)


# --------------------------------------------------------------------------
# 16.3 Average marker expression by cluster
# --------------------------------------------------------------------------

cat("Calculating average marker expression by cluster...\n")

average_expression <- AverageExpression(
  merfish_data,
  assays = "RNA",
  features = curated_markers,
  group.by = "seurat_clusters",
  slot = "data"
)

average_expression_matrix <- as.data.frame(
  average_expression$RNA
)

average_expression_matrix$gene <- rownames(
  average_expression_matrix
)


# Convert matrix to long format for plotting
average_expression_long <- average_expression_matrix %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "cluster",
    values_to = "expression"
  )


# Ensure clusters are ordered numerically
average_expression_long$cluster <- factor(
  average_expression_long$cluster,
  levels = paste0(
    "g",
    sort(
      as.numeric(
        gsub(
          "[^0-9]",
          "",
          unique(average_expression_long$cluster)
        )
      )
    )
  )
)


# In case Seurat returns numeric cluster names rather than g-prefixed names
if (!all(grepl("^g", as.character(average_expression_long$cluster)))) {

  average_expression_long$cluster <- factor(
    average_expression_long$cluster,
    levels = sort(
      unique(
        as.character(
          average_expression_long$cluster
        )
      )
    )
  )

}


# --------------------------------------------------------------------------
# 16.4 Average-expression heatmap
# --------------------------------------------------------------------------

cat("Generating average-expression marker heatmap...\n")

average_marker_heatmap <- ggplot(
  average_expression_long,
  aes(
    x = cluster,
    y = gene,
    fill = expression
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4C78A8",
    mid = "#F7F7D5",
    high = "#D73027",
    midpoint = 1
  ) +
  labs(
    title = "Average Marker Gene Expression Across Seurat Clusters",
    x = NULL,
    y = NULL,
    fill = "Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 10
    ),
    panel.grid = element_blank()
  )

ggsave(
  filename = "figures/marker_heatmap_average_expression.png",
  plot = average_marker_heatmap,
  width = 14,
  height = 10,
  dpi = 300
)

cat(
  "Marker heatmaps generated successfully.\n\n"
)


# ============================================================================
# 17. Manual Cluster Annotation
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 11: Cluster Annotation -----------------------------------------\n")
cat("=========================================================================\n")

cat("Applying marker-based manual cluster annotations...\n")


# Manual annotations were assigned after inspection of cluster marker
# profiles. Where the available marker profile did not support a precise
# subtype assignment, a broader or provisional identity was retained.

cluster_annotations <- data.frame(
  cluster = 0:19,

  cell_type = c(
    "Excitatory neuron (glutamatergic)",
    "GABAergic / inhibitory neuron",
    "Oligodendrocyte-lineage cell",
    "Excitatory neuron (glutamatergic)",
    "Astrocyte",
    "Endothelial cell",
    "Excitatory neuron (glutamatergic)",
    "GABAergic / inhibitory neuron",
    "Oligodendrocyte precursor cell (OPC)",
    "Striatal-like inhibitory neuron",
    "Microglia",
    "Vascular mural cell / pericyte",
    "Astrocyte",
    "Vascular smooth muscle cell",
    "Cholinergic neuron",
    "Oligodendrocyte-lineage / vascular-associated mixed cluster",
    "GABAergic / inhibitory neuron",
    "Excitatory neuron (glutamatergic)",
    "Neuronal subtype",
    "Oligodendrocyte-lineage cell"
  ),

  marker_genes = c(
    "Slc17a7, Grm2, Grin2b, Chrm1",
    "Slc32a1, Gpr101, Dlk1, Chrm5, Galr1",
    "S1pr5, Erbb3, Gjc3, Olig1",
    "Slc17a7, Baiap2, Chrm1, Celsr2",
    "Aqp4, Aldh1l1, Gfap, Sox9",
    "Cldn5, Flt1, Kdr, Tek",
    "Slc17a6, Grm4, Grm1, Ret",
    "Gad1, Slc32a1, Erbb4, Cnr1",
    "Gpr17, Pdgfra, Olig1, Sox8",
    "Adora2a, Drd2, Gpr88, Drd1",
    "Cx3cr1, P2ry12, Csf1r, C1qa",
    "Abcc9, Kcnj8, Pdgfrb, Myh11",
    "Aqp4, Gfap, Sox9, Aldh1l1",
    "Lmod1, Myh11, Hcar1, Pln",
    "Chat, Gpr151, Htr5b, Tacr1",
    "S1pr5, Erbb3, Gjc3, Olig1",
    "Cckar, Gad1, Slc32a1, Pln, Ret",
    "Slc17a7, Htr5b, Grin2b, Chrm1, Grm5",
    "Crhr2, F2rl3, Htr2c, Fzd4, Mc3r",
    "Gjc3, Olig1, S1pr5, Erbb3, Sox8"
  ),

  notes = c(
    "Strong Slc17a7 expression; glutamatergic neuronal profile",
    "Neuronal cluster with Slc32a1 expression and a distinctive receptor/neuropeptide-rich profile; precise subtype remains provisional",
    "Strong oligodendrocyte-lineage signature with Olig1, Sox8, Erbb3, S1pr5 and Gjc3",
    "Strong Slc17a7 expression with glutamatergic neuronal markers",
    "Strong astrocytic signature including Aqp4, Aldh1l1, Gfap and Sox9",
    "Strong endothelial signature including Cldn5, Flt1, Kdr, Tek, Emcn and Tie1",
    "Strong Slc17a6 expression with glutamatergic receptor and neuronal markers",
    "Strong Gad1 and Slc32a1 expression with Erbb4 and Cnr1, supporting an inhibitory neuronal identity",
    "Strong Gpr17, Pdgfra, Olig1 and Sox8 signature consistent with oligodendrocyte precursor cells",
    "Dopamine-receptor-rich inhibitory neuronal profile with Adora2a, Drd2, Gpr88 and Drd1; D1/D2 subtype remains unresolved",
    "Strong microglial signature including Cx3cr1, P2ry12, Csf1r, C1qa and C1qb",
    "Strong Abcc9/Kcnj8/Pdgfrb profile supports a vascular mural cell/pericyte identity; Myh11 suggests smooth-muscle-like characteristics",
    "Strong astrocytic profile with Aqp4, Gfap, Sox9 and Aldh1l1",
    "Strong contractile smooth-muscle profile dominated by Lmod1 and Myh11; vascular-associated genes are also present",
    "Strong Chat expression supports a cholinergic neuronal identity; Gpr151 and Htr5b suggest a specialized neuronal subtype",
    "Oligodendrocyte-lineage markers including S1pr5, Erbb3, Gjc3 and Olig1 are present alongside vascular-associated markers; identity remains provisional",
    "Strong Gad1, Slc32a1 and Erbb4 expression supports a GABAergic/inhibitory neuronal identity; Cckar and Pln indicate a specialized neuronal subtype, but precise subtype remains unresolved",
    "Strong Slc17a7 expression supports a glutamatergic neuronal identity; Htr5b, Gpr161, Drd5 and serotonin/dopamine receptor enrichment suggest a specialized neuronal subtype, but precise subtype remains unresolved",
    "Strong receptor-rich signature dominated by Crhr2, F2rl3 and Htr2c; precise neuronal subtype remains unresolved from the current marker profile",
    "Strong oligodendrocyte-lineage signature with Gjc3, Olig1, S1pr5, Erbb3 and Sox8; Myh11 suggests some vascular-associated signal, so precise subtype remains provisional"
  ),

  stringsAsFactors = FALSE
)


# Save annotation table
write.csv(
  cluster_annotations,
  file = "results/cluster_annotations_manual.csv",
  row.names = FALSE
)


# Save provisional annotation table separately
write.csv(
  cluster_annotations,
  file = "results/preliminary_cluster_annotations.csv",
  row.names = FALSE
)


# Add annotations to Seurat metadata
annotation_lookup <- setNames(
  cluster_annotations$cell_type,
  cluster_annotations$cluster
)

merfish_data$cell_type <- annotation_lookup[
  as.character(merfish_data$seurat_clusters)
]


# Create annotated UMAP
annotated_umap <- DimPlot(
  merfish_data,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("MERFISH UMAP - Annotated Cell Populations") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "figures/umap_annotated.png",
  plot = annotated_umap,
  width = 9,
  height = 7,
  dpi = 300
)

cat(
  "Manual annotations added to Seurat metadata.\n\n"
)


# ============================================================================
# 18. Spatial Coordinate Audit
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 12: Spatial Coordinate Audit ------------------------------------\n")
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


# ============================================================================
# 19. Save Processed and Annotated Seurat Objects
# ============================================================================

cat("=========================================================================\n")
cat("--- Step 13: Save Seurat Objects ----------------------------------------\n")
cat("=========================================================================\n")


# Save processed object before manual annotation
saveRDS(
  merfish_data,
  file = "results/merfish_processed.rds"
)


# Save final object containing cluster annotations
saveRDS(
  merfish_data,
  file = "results/merfish_seurat_annotated.rds"
)

cat(
  "Processed Seurat object saved to:\n",
  "results/merfish_processed.rds\n\n",
  sep = ""
)

cat(
  "Annotated Seurat object saved to:\n",
  "results/merfish_seurat_annotated.rds\n\n",
  sep = ""
)


# ============================================================================
# 20. Reproducibility Information
# ============================================================================

cat("=========================================================================\n")
cat("--- Reproducibility Information -----------------------------------------\n")
cat("=========================================================================\n")

capture.output(
  sessionInfo(),
  file = "results/sessionInfo.txt"
)

cat(
  "Session information saved to results/sessionInfo.txt\n\n"
)


# ============================================================================
# 21. Analysis Complete
# ============================================================================

cat("=========================================================================\n")
cat("--- Analysis Complete ---------------------------------------------------\n")
cat("=========================================================================\n\n")

cat("Generated outputs:\n\n")

cat("Figures:\n")
cat("  figures/qc_nFeature_RNA.png\n")
cat("  figures/pca_clusters.png\n")
cat("  figures/umap_clusters.png\n")
cat("  figures/spatial_clusters.png\n")
cat("  figures/marker_heatmap.png\n")
cat("  figures/marker_heatmap_average_expression.png\n")
cat("  figures/umap_annotated.png\n\n")

cat("Results:\n")
cat("  results/marker_genes_all_clusters.csv\n")
cat("  results/cluster_markers.csv\n")
cat("  results/top3_markers_per_cluster.csv\n")
cat("  results/top10_markers_per_cluster.csv\n")
cat("  results/curated_heatmap_markers.csv\n")
cat("  results/cluster_summary.csv\n")
cat("  results/cluster_annotations_manual.csv\n")
cat("  results/preliminary_cluster_annotations.csv\n")
cat("  results/qc_summary.csv\n")
cat("  results/spatial_coordinates.csv\n")
cat("  results/merfish_processed.rds\n")
cat("  results/merfish_seurat_annotated.rds\n")
cat("  results/sessionInfo.txt\n\n")

cat("MERFISH spatial transcriptomics analysis completed successfully.\n")

# ==============================================================================
# Single-cell RNA-seq analysis of primary and cultured human AT2 cells
# Samples: SB001 (P0), SB002 (P5), SB003 (P12)
# Clustering resolution: 0.1
# ==============================================================================

# -----------------------------
# 0. Load packages
# -----------------------------
library(Seurat)
library(sctransform)
library(matrixStats)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)

# Set the working directory to the analysis directory before running this script.


# -----------------------------
# 1. Read Cell Ranger output
# -----------------------------
# Update the paths below according to the local directory structure.

counts_001 <- Read10X(
  data.dir = "SB001_cellranger_count_outs/filtered_feature_bc_matrix/"
)
obj_001 <- CreateSeuratObject(
  counts = counts_001,
  project = "SB001",
  min.cells = 3,
  min.features = 200
)
obj_001$sample <- "SB001"

counts_002 <- Read10X(
  data.dir = "SB002_cellranger_count_outs/filtered_feature_bc_matrix/"
)
obj_002 <- CreateSeuratObject(
  counts = counts_002,
  project = "SB002",
  min.cells = 3,
  min.features = 200
)
obj_002$sample <- "SB002"

counts_003 <- Read10X(
  data.dir = "SB003_cellranger_count_outs/filtered_feature_bc_matrix/"
)
obj_003 <- CreateSeuratObject(
  counts = counts_003,
  project = "SB003",
  min.cells = 3,
  min.features = 200
)
obj_003$sample <- "SB003"

# -----------------------------
# 2. Merge samples
# -----------------------------
merged_obj <- merge(
  x = obj_001,
  y = list(obj_002, obj_003),
  add.cell.ids = c("SB001", "SB002", "SB003"),
  project = "AT2_culture"
)

# -----------------------------
# 3. Calculate mitochondrial percentage and visualize QC metrics
# -----------------------------
merged_obj[["percent.mt"]] <- PercentageFeatureSet(
  merged_obj,
  pattern = "^MT-"
)

p_qc_vln <- VlnPlot(
  merged_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "sample",
  ncol = 3
)

ggsave(
  filename = "QC_violin_plot.png",
  plot = p_qc_vln,
  width = 12,
  height = 5,
  dpi = 300
)

p_scatter_count_feature <- FeatureScatter(
  merged_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  filename = "QC_scatter_nCount_vs_nFeature.png",
  plot = p_scatter_count_feature,
  width = 6,
  height = 5,
  dpi = 300
)

p_scatter_count_mt <- FeatureScatter(
  merged_obj,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

ggsave(
  filename = "QC_scatter_nCount_vs_percentMT.png",
  plot = p_scatter_count_mt,
  width = 6,
  height = 5,
  dpi = 300
)

# -----------------------------
# 4. QC summary
# -----------------------------
qc_summary <- merged_obj@meta.data %>%
  group_by(sample) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA),
    mean_nFeature_RNA = mean(nFeature_RNA),
    median_nCount_RNA = median(nCount_RNA),
    mean_nCount_RNA = mean(nCount_RNA),
    median_percent_mt = median(percent.mt),
    mean_percent_mt = mean(percent.mt),
    pct_high_mt_10 = mean(percent.mt > 10) * 100,
    pct_high_mt_15 = mean(percent.mt > 15) * 100,
    pct_high_mt_20 = mean(percent.mt > 20) * 100,
    pct_low_features_200 = mean(nFeature_RNA < 200) * 100,
    pct_low_features_500 = mean(nFeature_RNA < 500) * 100,
    pct_high_features_6000 = mean(nFeature_RNA > 6000) * 100,
    pct_high_features_7500 = mean(nFeature_RNA > 7500) * 100
  )

write.xlsx(
  qc_summary,
  file = "QC_summary_SB001_SB003.xlsx",
  rowNames = FALSE
)

# -----------------------------
# 5. Quality control filtering
# -----------------------------
# Cells were retained based on the following criteria:
# - nFeature_RNA > 500
# - nFeature_RNA < 8,000
# - percent.mt < 15%

merged_obj_qc <- subset(
  merged_obj,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 8000 &
    percent.mt < 15
)

# Summarize cell retention after QC filtering.
before_qc <- table(merged_obj$sample)
after_qc <- table(merged_obj_qc$sample)

qc_retention <- data.frame(
  sample = names(before_qc),
  before_qc = as.numeric(before_qc),
  after_qc = as.numeric(after_qc),
  retained_percent = as.numeric(after_qc) / as.numeric(before_qc) * 100
)

write.xlsx(
  qc_retention,
  file = "QC_retention_SB001_SB003.xlsx",
  rowNames = FALSE
)

# -----------------------------
# 6. SCTransform normalization
# -----------------------------
# Normalize UMI counts using SCTransform and regress out the percentage of
# mitochondrial transcripts.
plan(sequential)
options(future.globals.maxSize = 5000 * 1024^2)

merged_obj_qc <- SCTransform(
  merged_obj_qc,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

# -----------------------------
# 7. Principal component analysis
# -----------------------------
merged_obj_qc <- RunPCA(
  merged_obj_qc,
  verbose = FALSE
)

p_pca_sample <- DimPlot(
  merged_obj_qc,
  reduction = "pca",
  group.by = "sample"
)

ggsave(
  filename = "PCA_by_sample.png",
  plot = p_pca_sample,
  width = 6,
  height = 5,
  dpi = 300
)

# Inspect the elbow plot to guide selection of the number of PCs.
p_elbow <- ElbowPlot(
  merged_obj_qc,
  ndims = 50
)

ggsave(
  filename = "PCA_elbow_plot.png",
  plot = p_elbow,
  width = 6,
  height = 5,
  dpi = 300
)

DefaultAssay(merged_obj_qc) <- "SCT"
dims_use <- 1:20

# -----------------------------
# 8. UMAP dimensional reduction
# -----------------------------
merged_obj_qc <- RunUMAP(
  merged_obj_qc,
  dims = dims_use
)

# Assign passage labels.
merged_obj_qc$condition <- recode(
  merged_obj_qc$sample,
  "SB001" = "P0",
  "SB002" = "P5",
  "SB003" = "P12"
)
merged_obj_qc$condition <- factor(
  merged_obj_qc$condition,
  levels = c("P0", "P5", "P12")
)

p_umap_condition <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "condition"
) +
  coord_fixed() +
  theme_classic() +
  ggtitle("UMAP by passage")

# -----------------------------
# 9. Graph-based clustering
# -----------------------------
merged_obj_qc <- FindNeighbors(
  merged_obj_qc,
  dims = dims_use
)

merged_obj_qc <- FindClusters(
  merged_obj_qc,
  resolution = 0.1
)

p_umap_cluster <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP by cluster")

p_umap_combined <- p_umap_condition + p_umap_cluster

ggsave(
  filename = "UMAP_passage_and_cluster.png",
  plot = p_umap_combined,
  width = 12,
  height = 5,
  dpi = 300
)

# -----------------------------
# 10. Cluster composition by sample
# -----------------------------
cluster_sample_count <- as.data.frame(
  table(
    Sample = merged_obj_qc$sample,
    Cluster = merged_obj_qc$seurat_clusters
  )
)
colnames(cluster_sample_count) <- c("Sample", "Cluster", "Cell_count")

cluster_sample_percent <- as.data.frame(
  prop.table(
    table(
      Sample = merged_obj_qc$sample,
      Cluster = merged_obj_qc$seurat_clusters
    ),
    margin = 1
  ) * 100
)
colnames(cluster_sample_percent) <- c("Sample", "Cluster", "Percent")

write.xlsx(
  list(
    "Cell_count" = cluster_sample_count,
    "Percent_by_sample" = cluster_sample_percent
  ),
  file = "cluster_sample_composition.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)

# Use the in-memory object directly rather than reloading the exported Excel file.
df_cluster_percent <- cluster_sample_percent %>%
  mutate(
    Sample = recode(
      Sample,
      "SB001" = "P0",
      "SB002" = "P5",
      "SB003" = "P12"
    ),
    Sample = factor(Sample, levels = c("P0", "P5", "P12")),
    Cluster = factor(Cluster)
  )

p_stacked <- ggplot(
  df_cluster_percent,
  aes(x = Sample, y = Percent, fill = Cluster)
) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(
    x = "Sample",
    y = "Percentage of cells (%)",
    fill = "Cluster"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = "P0_P5_P12_cluster_percent_stacked_barplot.png",
  plot = p_stacked,
  width = 6,
  height = 5,
  dpi = 300
)


# -----------------------------
# 11. Feature plots for known markers
# -----------------------------
DefaultAssay(merged_obj_qc) <- "SCT"

# Common plotting theme used for marker FeaturePlots.
feature_plot_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 15, face = "italic", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(1, 1, 10, 1),
    legend.margin = margin(t = 0, r = 10, b = 0, l = -5)
  )

feature_plot_guides <- guides(
  color = guide_colorbar(
    barheight = unit(0.5, "cm"),
    barwidth = unit(0.25, "cm"),
    title.position = "top"
  )
)

# AT2 markers
p_AT2 <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "NKX2-1", "SFTPC", "SFTPA1", "SFTPB",
    "SFTPD", "SLC34A2", "ABCA3", "NAPSA"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_AT2_markers.png",
  plot = p_AT2,
  width = 14,
  height = 7,
  dpi = 300
)

# Progenitor-like AT2 markers
p_progenitorAT2 <- FeaturePlot(
  merged_obj_qc,
  features = c("TM4SF1", "AXIN2", "ETV5", "FGFR2", "WIF1", "NKD1"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_progenitorAT2_markers.png",
  plot = p_progenitorAT2,
  width = 14,
  height = 7,
  dpi = 300
)

# AT1-like markers
p_AT1 <- FeaturePlot(
  merged_obj_qc,
  features = c("AGER", "PDPN", "AQP5", "CAV1", "CLIC5"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_AT1_like_markers.png",
  plot = p_AT1,
  width = 14,
  height = 7,
  dpi = 300
)

# Basal cell markers
p_Basal <- FeaturePlot(
  merged_obj_qc,
  features = c("KRT5", "KRT17"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_Basal_cell_markers.png",
  plot = p_Basal,
  width = 14,
  height = 7,
  dpi = 300
)

# Proliferation markers
p_proliferate <- FeaturePlot(
  merged_obj_qc,
  features = c("MKI67", "TOP2A", "UBE2C", "CENPF"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_Proliferation_markers.png",
  plot = p_proliferate,
  width = 14,
  height = 7,
  dpi = 300
)

# Transitional-state markers
p_Trans <- FeaturePlot(
  merged_obj_qc,
  features = c("KRT8", "KRT7", "CLDN4", "KRT19", "KRT17", "MMP7", "LGALS3"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_Transitional_markers.png",
  plot = p_Trans,
  width = 14,
  height = 7,
  dpi = 300
)

# Stress/injury-associated markers
p_Stress <- FeaturePlot(
  merged_obj_qc,
  features = c("ATF3", "FOS", "JUN", "HSPA1A", "DDIT3", "CDKN1A"),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  feature_plot_theme &
  feature_plot_guides

ggsave(
  filename = "FeaturePlot_Stress_markers.png",
  plot = p_Stress,
  width = 14,
  height = 7,
  dpi = 300
)

# -----------------------------
# 12. AT2 gene signature module score
# -----------------------------
AT2_genes <- c(
  "SFTPC", "SFTPB", "SFTPA1", "SFTPA2",
  "ABCA3", "NAPSA", "LAMP3", "SLC34A2"
)

signature_list <- list(
  AT2 = AT2_genes
)

# Retain only genes present in the Seurat object.
signature_list_use <- lapply(signature_list, function(genes) {
  genes[genes %in% rownames(merged_obj_qc)]
})

DefaultAssay(merged_obj_qc) <- "SCT"

merged_obj_qc <- AddModuleScore(
  object = merged_obj_qc,
  features = signature_list_use,
  name = names(signature_list_use)
)

AT2signature <- FeaturePlot(
  merged_obj_qc,
  features = "AT21",
  reduction = "umap"
)

ggsave(
  filename = "FeaturePlot_AT2_signature.png",
  plot = AT2signature,
  width = 14,
  height = 7,
  dpi = 300
)

# -----------------------------
# 13. Identify cluster marker genes
# -----------------------------
DefaultAssay(merged_obj_qc) <- "SCT"
Idents(merged_obj_qc) <- "seurat_clusters"

# Prepare SCT models for differential expression testing.
merged_obj_qc <- PrepSCTFindMarkers(merged_obj_qc)

cluster_markers <- FindAllMarkers(
  merged_obj_qc,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.05
)

# Define a more stringent subset of cluster markers.
# pct.1: fraction of cells expressing the gene in the target cluster.
# pct.2: fraction of cells expressing the gene in all other cells.
cluster_markers_strong <- cluster_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > 0.10
  )

write.xlsx(
  cluster_markers,
  file = "cluster_markers_FindAllMarkers.xlsx",
  rowNames = FALSE
)

write.xlsx(
  cluster_markers_strong,
  file = "cluster_markers_strong.xlsx",
  rowNames = FALSE
)

# -----------------------------
# 14. Select top markers per cluster
# -----------------------------
top20_markers <- cluster_markers_strong %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()

top50_markers <- cluster_markers_strong %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50, with_ties = FALSE) %>%
  ungroup()

# -----------------------------
# 15. Heatmaps of top cluster markers
# -----------------------------
DefaultAssay(merged_obj_qc) <- "SCT"
Idents(merged_obj_qc) <- "seurat_clusters"

# Downsample to a maximum of 300 cells per cluster for visualization.
set.seed(123)

cells_use <- unlist(
  lapply(levels(Idents(merged_obj_qc)), function(clust) {
    cells <- WhichCells(merged_obj_qc, idents = clust)
    sample(cells, size = min(300, length(cells)))
  })
)

p_heat_top50 <- DoHeatmap(
  merged_obj_qc,
  features = unique(top50_markers$gene),
  cells = cells_use,
  group.by = "seurat_clusters",
  size = 3
) +
  NoLegend()

ggsave(
  filename = "cluster_top50_marker_heatmap.pdf",
  plot = p_heat_top50,
  width = 12,
  height = 20
)

p_heat_top20 <- DoHeatmap(
  merged_obj_qc,
  features = unique(top20_markers$gene),
  cells = cells_use,
  group.by = "seurat_clusters",
  size = 3
) +
  NoLegend()

ggsave(
  filename = "cluster_top20_marker_heatmap_res.pdf",
  plot = p_heat_top20,
  width = 12,
  height = 14
)

# -----------------------------
# 16. Dot plots for marker genes
# -----------------------------
top20_genes <- unique(top20_markers$gene)

p_dot <- DotPlot(
  merged_obj_qc,
  features = top20_genes,
  group.by = "seurat_clusters"
) +
  scale_color_gradient(
    low = "grey90",
    high = "purple"
  ) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  filename = "cluster_top20_marker_dotplot_res.pdf",
  plot = p_dot,
  width = 16,
  height = 8
)

# Selected known marker genes.
genes_to_plot <- c(
  "SFTPC", "SFTPA1", "SFTPB", "SFTPD", "SLC34A2", "NAPSA", "NKX2-1",
  "AGER", "PDPN", "CLIC5", "CAV1",
  "KRT5", "KRT17",
  "SOX9", "COL1A1", "CDKN2B", "EPCAM",
  "MKI67", "TOP2A"
)

p_dot_selected <- DotPlot(
  merged_obj_qc,
  features = genes_to_plot,
  group.by = "seurat_clusters",
  dot.scale = 8
) +
  scale_color_gradient(
    low = "grey90",
    high = "purple"
  ) +
  RotatedAxis()

ggsave(
  filename = "cluster_marker_selected_dotplot.pdf",
  plot = p_dot_selected,
  width = 16,
  height = 6
)

p_dot_selected_by_sample <- DotPlot(
  merged_obj_qc,
  features = genes_to_plot,
  group.by = "sample",
  dot.scale = 8
) +
  scale_color_gradient(
    low = "grey90",
    high = "purple"
  ) +
  RotatedAxis()

ggsave(
  filename = "cluster_marker_selected_dotplot_by_sample.pdf",
  plot = p_dot_selected_by_sample,
  width = 16,
  height = 3
)

# -----------------------------
# 17. Violin plots by marker group
# -----------------------------
DefaultAssay(merged_obj_qc) <- "SCT"

AT2_markers <- c(
  "SFTPC", "SFTPA1", "SFTPB", "SFTPA2",
  "SFTPD", "SLC34A2", "ABCA3", "NAPSA"
)

AT1_markers <- c(
  "AGER", "PDPN", "AQP5", "CAV1",
  "CAV2", "CLIC5", "HOPX"
)

Basal_markers <- c(
  "KRT5", "TP63"
)

Transition_markers <- c(
  "CLDN4", "KRT7", "KRT8", "KRT17", "KRT19"
)

# Retain only markers present in the dataset.
AT2_present <- intersect(AT2_markers, rownames(merged_obj_qc))
AT1_present <- intersect(AT1_markers, rownames(merged_obj_qc))
Basal_present <- intersect(Basal_markers, rownames(merged_obj_qc))
Transition_present <- intersect(Transition_markers, rownames(merged_obj_qc))

violin_theme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  plot.title = element_text(
    face = "italic",
    size = 20
  )
)

# AT2 markers
p_vln_AT2 <- VlnPlot(
  merged_obj_qc,
  features = AT2_present,
  group.by = "condition",
  pt.size = 0,
  ncol = 4
) &
  violin_theme

ggsave(
  filename = "Violin_AT2_markers.pdf",
  plot = p_vln_AT2,
  width = 10,
  height = 6
)
ggsave(
  filename = "Violin_AT2_markers.png",
  plot = p_vln_AT2,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# AT1 markers
p_vln_AT1 <- VlnPlot(
  merged_obj_qc,
  features = AT1_present,
  group.by = "condition",
  pt.size = 0,
  ncol = 4
) &
  violin_theme

ggsave(
  filename = "Violin_AT1_markers.pdf",
  plot = p_vln_AT1,
  width = 10,
  height = 6
)
ggsave(
  filename = "Violin_AT1_markers.png",
  plot = p_vln_AT1,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Basal cell markers
p_vln_Basal <- VlnPlot(
  merged_obj_qc,
  features = Basal_present,
  group.by = "condition",
  pt.size = 0,
  ncol = 4
) &
  violin_theme

ggsave(
  filename = "Violin_Basal_markers.pdf",
  plot = p_vln_Basal,
  width = 10,
  height = 3
)
ggsave(
  filename = "Violin_Basal_markers.png",
  plot = p_vln_Basal,
  width = 10,
  height = 3,
  dpi = 300,
  bg = "white"
)

# Transitional-state markers
p_vln_Transition <- VlnPlot(
  merged_obj_qc,
  features = Transition_present,
  group.by = "condition",
  pt.size = 0,
  ncol = 4
) &
  violin_theme

ggsave(
  filename = "Violin_Transition_markers.pdf",
  plot = p_vln_Transition,
  width = 10,
  height = 3
)
ggsave(
  filename = "Violin_Transition_markers.png",
  plot = p_vln_Transition,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

# -----------------------------
# 18. Save processed Seurat object
# -----------------------------
saveRDS(
  merged_obj_qc,
  file = "SB001_SB003_processed.rds"
)


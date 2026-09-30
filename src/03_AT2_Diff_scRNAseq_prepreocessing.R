# ==============================================================================
# Single-cell RNA-seq analysis of differentiated human AT2 cells
# Conditions: DMSO, LATSi, MAPKi, LATSi+MAPKi
# Samples: SB021, SB022, SB023, SB024
# Final clustering resolution: 0.1
# ==============================================================================

library(Seurat)
library(sctransform)
library(matrixStats)
library(future)
library(ggplot2)
library(patchwork)
library(dplyr)
library(openxlsx)
library(pheatmap)

# =========================
# 1. Read Cell Ranger output
# =========================
# Use the filtered feature-barcode matrix from each Cell Ranger sample output.
# Update paths to match the local directory structure.

# Set the working directory to the differentiation analysis directory before running.
# setwd("path/to/Diff")

# SB021
counts_021 <- Read10X(
  data.dir = "SB2124_cellranger_multi_outs/per_sample_outs/SB021/count/sample_filtered_feature_bc_matrix/"
)

obj_021 <- CreateSeuratObject(
  counts = counts_021,
  project = "SB021",
  min.cells = 3,
  min.features = 200
)

obj_021$sample <- "SB021"


# SB022
counts_022 <- Read10X(
  data.dir = "SB2124_cellranger_multi_outs/per_sample_outs/SB022/count/sample_filtered_feature_bc_matrix/"
)

obj_022 <- CreateSeuratObject(
  counts = counts_022,
  project = "SB022",
  min.cells = 3,
  min.features = 200
)

obj_022$sample <- "SB022"


# SB023
counts_023 <- Read10X(
  data.dir = "SB2124_cellranger_multi_outs/per_sample_outs/SB023/count/sample_filtered_feature_bc_matrix/"
)

obj_023 <- CreateSeuratObject(
  counts = counts_023,
  project = "SB023",
  min.cells = 3,
  min.features = 200
)

obj_023$sample <- "SB023"

# SB024
counts_024 <- Read10X(
  data.dir = "SB2124_cellranger_multi_outs/per_sample_outs/SB024/count/sample_filtered_feature_bc_matrix/"
)

obj_024 <- CreateSeuratObject(
  counts = counts_024,
  project = "SB024",
  min.cells = 3,
  min.features = 200
)

obj_024$sample <- "SB024"

# =========================
# 2. Merge samples
# =========================

merged_obj <- merge(
  x = obj_021,
  y = list(obj_022, obj_023, obj_024),
  add.cell.ids = c("SB021", "SB022", "SB023", "SB024"),
  project = "AT2_Diff"
)


# =========================
# 3. Calculate mitochondrial percentage
# =========================

merged_obj[["percent.mt"]] <- PercentageFeatureSet(
  merged_obj,
  pattern = "^MT-"
)

# QC violin plot
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


# Feature scatter plot

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

# nCount_RNA vs mitochondrial percentage
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

# =========================
# 4. QC summary table
# =========================

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
  file = "QC_summary_SB021_SB024.xlsx",
  rowNames = FALSE
)

# =========================
# 5. QC filtering
# =========================
# Filtering criteria:
# nFeature_RNA > 500
# nFeature_RNA < 8000
# percent.mt < 15

merged_obj_qc <- subset(
  merged_obj,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 8000 &
    percent.mt < 15
)

# Compare cell numbers before and after QC filtering.
before_qc <- table(merged_obj$sample)
after_qc <- table(merged_obj_qc$sample)

qc_retention <- data.frame(
  sample = names(before_qc),
  before_qc = as.numeric(before_qc),
  after_qc = as.numeric(after_qc),
  retained_percent = as.numeric(after_qc) / as.numeric(before_qc) * 100
)


# =========================
# 6. SCTransform normalization
# Normalize UMI counts and regress out mitochondrial percentage.
# =========================

plan(sequential)
options(future.globals.maxSize = 5000 * 1024^2)

DefaultAssay(merged_obj_qc) <- "RNA"
merged_obj_qc <- JoinLayers(merged_obj_qc)

merged_obj_qc <- SCTransform(
  merged_obj_qc,
  vars.to.regress = "percent.mt",
  verbose = TRUE
)


# =========================
# 7. PCA
# =========================

merged_obj_qc <- RunPCA(
  merged_obj_qc,
  verbose = FALSE
)

# Inspect sample distribution in PCA space.
DimPlot(
  merged_obj_qc,
  reduction = "pca",
  group.by = "sample"
)

# Inspect the elbow plot to select the number of PCs.
ElbowPlot(
  merged_obj_qc,
  ndims = 50
)

# =========================
# 8. Set PCs to use
# =========================
DefaultAssay(merged_obj_qc) <- "SCT"

dims_use <- 1:20

# =========================
# 9. UMAP
# =========================

merged_obj_qc <- RunUMAP(
  merged_obj_qc,
  dims = dims_use
)

# Inspect UMAP by sample.
DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "sample"
)

# Optional UMAP parameter adjustment
# merged_obj_qc <- RunUMAP(
#   merged_obj_qc,
#   dims = dims_use,
#   min.dist = 0.1,
#   spread = 1
# )
# 
# DimPlot(
#   merged_obj_qc,
#   reduction = "umap",
#   group.by = "sample"
# )

# =========================
# 10. Clustering
# =========================
# Add experimental condition metadata.
merged_obj_qc$condition <- merged_obj_qc$sample

merged_obj_qc$condition[merged_obj_qc$sample == "SB021"] <- "DMSO"
merged_obj_qc$condition[merged_obj_qc$sample == "SB022"] <- "LATSi"
merged_obj_qc$condition[merged_obj_qc$sample == "SB023"] <- "MAPKi"
merged_obj_qc$condition[merged_obj_qc$sample == "SB024"] <- "LATSi+MAPKi"

condition_order <- c(
  "DMSO",
  "LATSi",
  "MAPKi",
  "LATSi+MAPKi"
)

merged_obj_qc$condition <- factor(
  merged_obj_qc$condition,
  levels = condition_order
)


# UMAP by condition
p_umap_condition <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "condition"
) +
  coord_fixed() +
  theme_classic()


ggsave(
  filename = "UMAP_by_condition.png",
  plot = p_umap_condition,
  width = 6,
  height = 5,
  dpi = 300
)


merged_obj_qc <- FindNeighbors(
  merged_obj_qc,
  dims = dims_use
)

merged_obj_qc <- FindClusters(
  merged_obj_qc,
  resolution = 0.1
)

DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)

##
library(ggplot2)
library(R.utils)
library(Matrix)

p_umap_cluster <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP by cluster")

# Cluster cell counts by sample
# table(merged_obj_qc$sample, merged_obj_qc$seurat_clusters)
# 
# Cluster proportions within each sample
# prop.table(
#   table(merged_obj_qc$sample, merged_obj_qc$seurat_clusters),
#   margin = 1
# ) * 100

# Cluster cell counts by sample

cluster_sample_count <- as.data.frame(
  table(
    Sample = merged_obj_qc$sample,
    Cluster = merged_obj_qc$seurat_clusters
  )
)

colnames(cluster_sample_count) <- c("Sample", "Cluster", "Cell_count")


# Cluster proportions within each sample

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


# Save cluster counts and proportions.
write.xlsx(
  list(
    "Cell_count" = cluster_sample_count,
    "Percent_by_sample" = cluster_sample_percent
  ),
  file = "cluster_sample_composition.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)


# =========================
# 11. UMAP by condition and cluster
# =========================


# ============================================================
# Common UMAP theme
# ============================================================

umap_theme <- theme_classic() +
  theme(
    aspect.ratio = 1,
    
    plot.title = element_text(
      size = 15,
      hjust = 0.5
    ),
    
    axis.title.x = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 4)
    ),
    
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      margin = margin(r = 4)
    ),
    
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.4,
      linetype = "solid"
    ),
    
    axis.line = element_blank(),
    
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 15)
  )


# ============================================================
# UMAP by condition
# ============================================================


p_umap_condition <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "condition"
) +
  scale_color_discrete(
    breaks = c("DMSO", "LATSi", "MAPKi", "LATSi+MAPKi")
  ) +
  ggtitle("UMAP by condition") +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  coord_fixed() +
  umap_theme


# ============================================================
# UMAP by cluster
# ============================================================
cluster_colors <- c(
  "0" = "#1F77B4",  # blue
  "1" = "#FF7F0E",  # orange
  "2" = "#2CA02C",  # green
  "3" = "#D62728",  # red
  "4" = "#9467BD",  # purple
  "5" = "#8C564B",  # brown
  "6" = "#E377C2",  # pink
  "7" = "#7F7F7F",  # gray
  "8" = "#BCBD22",  # olive
  "9" = "#17BECF"   # cyan
)

# cluster_colors <- c(
#   "0" = "#4E79A7",
#   "1" = "#F28E2B",
#   "2" = "#E15759",
#   "3" = "#76B7B2",
#   "4" = "#59A14F",
#   "5" = "#EDC948",
#   "6" = "#B07AA1",
#   "7" = "#FF9DA7",
#   "8" = "#9C755F",
#   "9" = "#BAB0AC"
# )
p_umap_cluster <- DimPlot(
  merged_obj_qc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) +
  scale_color_manual(
    values = cluster_colors
  ) +
  ggtitle("UMAP by cluster") +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  coord_fixed() +
  umap_theme


# ============================================================
# Combine UMAP panels
# ============================================================

p_umap_combined <- p_umap_condition + p_umap_cluster

p_umap_combined


# ============================================================
# Save combined UMAP
# ============================================================

ggsave(
  filename = "UMAP_condition_and_cluster.png",
  plot = p_umap_combined,
  width = 10,
  height = 5,
  dpi = 300
)
# =========================
# 12. FeaturePlot for known markers
# =========================

DefaultAssay(merged_obj_qc) <- "RNA"

# Merge RNA layers if needed
merged_obj_qc <- JoinLayers(
  merged_obj_qc,
  assay = "RNA"
)

# Create normalized "data" layer
merged_obj_qc <- NormalizeData(
  merged_obj_qc,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# Common theme for marker FeaturePlots
feature_theme <- theme_classic() +
  theme(
    aspect.ratio = 1,
    
    plot.title = element_text(
      size = 17,
      face = "italic",
      hjust = 0.5
    ),
    
    # axis titles
    axis.title.x = element_text(
      size = 10,
      face = "bold",
      margin = margin(t = 4)
    ),
    axis.title.y = element_text(
      size = 10,
      face = "bold",
      margin = margin(r = 4)
    ),
    
    # no tick labels
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    # solid border
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.4,
      linetype = "solid"
    ),
    
    # remove classic axis lines because panel.border is used
    axis.line = element_blank(),
    
    # legend
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(
      size = 7,
      margin = margin(t = 1)
    ),
    
    # remove legend tick marks
    legend.ticks = element_blank(),
    legend.axis.line = element_blank(),
    
    plot.margin = margin(3, 3, 3, 3),
    
    legend.margin = margin(
      t = -2,
      r = 0,
      b = 0,
      l = 0
    )
  )


feature_guide <- guides(
  color = guide_colorbar(
    direction = "horizontal",
    barwidth = grid::unit(2.5, "cm"),
    barheight = grid::unit(0.15, "cm"),
    title.position = "top",
    title.hjust = 0.5,
    ticks = FALSE
  )
)

# ============================================================
# Downstream gene markers
# ============================================================

p_ERK <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "IER3",
    "EGR1",
    "SPRY2",
    "SPRY4",
    "ETV5",
    "DUSP4",
    "DUSP5"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_ERK_markers.png",
  plot = p_ERK,
  width = 14,
  height = 8,
  dpi = 300
)

p_YAP <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "INHBA",
    "ANKRD1",
    "CCN1",
    "CCN2"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_YAP_markers.png",
  plot = p_YAP,
  width = 14,
  height = 8,
  dpi = 300
)

# ============================================================
# AT2 markers
# ============================================================

p_AT2 <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "NKX2-1",
    "SFTPC",
    "SFTPA1",
    "SFTPB",
    "SFTPD",
    "SLC34A2",
    "ABCA3",
    "NAPSA"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_AT2_markers.png",
  plot = p_AT2,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# Progenitor-like AT2 markers
# ============================================================

p_progenitorAT2 <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "TM4SF1",
    "AXIN2",
    "ETV5",
    "FGFR2",
    "WIF1",
    "NKD1"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_progenitorAT2_markers.png",
  plot = p_progenitorAT2,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# AT1-like markers
# ============================================================

p_AT1 <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "AGER",
    "PDPN",
    "AQP5",
    "CAV1",
    "CLIC5",
    "HOPX"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_AT1_like_markers.png",
  plot = p_AT1,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# Basal cell markers
# ============================================================

p_Basal <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "KRT5",
    "TP63"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_Basal_cell_markers.png",
  plot = p_Basal,
  width = 7,
  height = 4.5,
  dpi = 300
)


# ============================================================
# Proliferation markers
# ============================================================

p_proliferate <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_Proliferation_markers.png",
  plot = p_proliferate,
  width = 14,
  height = 4.5,
  dpi = 300
)


# ============================================================
# Transitional markers
# ============================================================

p_Trans <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "KRT8",
    "KRT7",
    "CLDN4",
    "KRT19",
    "KRT17",
    "MMP7",
    "LGALS3"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_Transitional_markers.png",
  plot = p_Trans,
  width = 14,
  height = 8,
  dpi = 300
)


# ============================================================
# Stress-injury markers
# ============================================================

p_Stress <- FeaturePlot(
  merged_obj_qc,
  features = c(
    "ATF3",
    "FOS",
    "JUN",
    "HSPA1A",
    "DDIT3",
    "CDKN1A"
  ),
  reduction = "umap",
  ncol = 4
) &
  coord_fixed() &
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) &
  feature_theme &
  feature_guide


ggsave(
  filename = "FeaturePlot_Stress_markers.png",
  plot = p_Stress,
  width = 14,
  height = 8,
  dpi = 300
)

# =========================
# 13. AT1 and AT2 gene signature module scores
# =========================
# =========================
# AT1 signature
# Burgess et al., Cell Stem Cell 2024
# AT1 vs all lung cells top 50 genes
# =========================

AT1_genes <- c(
  "AGER", "EMP2", "CAV1", "CEACAM6", "HOPX",
  "MYL9", "GPRC5A", "CLDN18", "ADIRF", "KRT7",
  "RTKN2", "LMO7", "SFTA2", "CYP4B1", "CLIC3",
  "C19orf33", "KRT19", "SPOCK2", "TSPAN13", "CAV2",
  "FXYD3", "TNNC1", "ANXA3", "SLC39A8", "AQP4",
  "TACSTD2", "SCEL", "RNASE1", "CST6", "SFTA1P",
  "KRT18", "FOLR1", "SUSD2", "CD55", "LIMCH1",
  "ATP5F1E", "ANKRD29", "MSLN", "EPCAM", "SELENOW",
  "IGFBP7", "KRT8", "VEGFA", "CLIC5", "ANKRD1",
  "SEMA3B", "GGTLC1", "BEX3", "CADM1", "PEBP4"
)

DefaultAssay(merged_obj_qc) <- "SCT"

# Retain genes present in the SCT assay.
AT1_genes_use <- AT1_genes[
  AT1_genes %in% rownames(merged_obj_qc)
]

merged_obj_qc <- AddModuleScore(
  object = merged_obj_qc,
  features = list(AT1_genes_use),
  name = "AT1_signature"
)

AT1signature <- FeaturePlot(
  merged_obj_qc,
  features = "AT1_signature1",
  reduction = "umap",
  pt.size = 0.2,
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster = FALSE
) +
  ggtitle("AT1 Signature") +
  coord_fixed() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  feature_theme +
  feature_guide


# PDF - publication/vector version
ggsave(
  "FeaturePlot_AT1_signature.pdf",
  AT1signature,
  width = 5,
  height = 5,
  device = cairo_pdf
)

# PNG - preview/image version
ggsave(
  "FeaturePlot_AT1_signature.png",
  AT1signature,
  width = 5,
  height = 5,
  dpi = 600
)

# =========================
# AT2 signature
# Burgess et al., Cell Stem Cell 2024
# AT2 vs all lung cells top 50 genes
# =========================

AT2_genes <- c(
  "SFTPA1",  "SFTPA2",  "NAPSA",  "SFTPD",  "SFTPB",
  "SFTA2",  "PGC",  "SFTPC",  "SLC34A2",  "PEBP4",
  "SLPI",  "MUC1",  "CXCL17",  "WIF1",  "HOPX",
  "LRRK2",  "LAMP3",  "S100A14",  "RNASE1",  "C16orf89",
  "FGGY",  "SDR16C5",  "ABCA3",  "SFTA3",  "C11orf96",
  "SELENBP1",  "PLA2G1B",  "MFSD2A",  "AK1",  "GKN2",
  "LPCAT1",  "HHIP",  "KRT18",  "SOD3",  "CEBPD",
  "PIGR",  "KRT8",  "SELENOP",  "MALL",  "C3",  "MSMO1",
  "CLDN18",  "SDC4",  "SLC39A8",  "CTSH",  "C4BPA",
  "RNF145",  "TFPI",  "RGS16",  "SLC22A31"
)

DefaultAssay(merged_obj_qc) <- "SCT"

# Retain genes present in the SCT assay.
AT2_genes_use <- AT2_genes[
  AT2_genes %in% rownames(merged_obj_qc)
]

merged_obj_qc <- AddModuleScore(
  object = merged_obj_qc,
  features = list(AT2_genes_use),
  name = "AT2_signature"
)

AT2signature <- FeaturePlot(
  merged_obj_qc,
  features = "AT2_signature1",
  reduction = "umap",
  pt.size = 0.2,
  order = TRUE,
  min.cutoff = "q10",
  max.cutoff = "q90",
  raster = FALSE
) +
  ggtitle("AT2 Signature") +
  coord_fixed() +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  feature_theme +
  feature_guide


# PDF - publication/vector version
ggsave(
  "FeaturePlot_AT2_signature.pdf",
  AT2signature,
  width = 5,
  height = 5,
  device = cairo_pdf
)

# PNG - preview/image version
ggsave(
  "FeaturePlot_AT2_signature.png",
  AT2signature,
  width = 5,
  height = 5,
  dpi = 600
)
# =========================
# 14. FindAllMarkers: cluster markers
# =========================

DefaultAssay(merged_obj_qc) <- "SCT"
Idents(merged_obj_qc) <- "seurat_clusters"

# Prepare the SCT assay for marker testing.
merged_obj_qc <- PrepSCTFindMarkers(merged_obj_qc)

cluster_markers <- FindAllMarkers(
  merged_obj_qc,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.05
)


# Filter for strong cluster markers.
# pct.1: fraction of cells expressing the gene in the target cluster.
# pct.2: fraction of cells expressing the gene outside the target cluster.

cluster_markers_strong <- cluster_markers %>%
  filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > 0.10
  )


# Save marker tables.
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

# =========================
# 15. Select top markers per cluster
# =========================

top20_markers <- cluster_markers_strong %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup()


# =========================
# 16. Heatmap for top marker genes
# =========================

DefaultAssay(merged_obj_qc) <- "SCT"
Idents(merged_obj_qc) <- "seurat_clusters"

# Sample up to 300 cells per cluster for the heatmap.
set.seed(123)

cells_use <- unlist(
  lapply(levels(Idents(merged_obj_qc)), function(clust) {
    cells <- WhichCells(merged_obj_qc, idents = clust)
    sample(cells, size = min(300, length(cells)))
  })
)

p_heatmap_top20 <- DoHeatmap(
  merged_obj_qc,
  features = unique(top20_markers$gene),
  cells = cells_use,
  group.by = "seurat_clusters",
  size = 3
) +
  NoLegend()

ggsave(
  filename = "cluster_top20_marker_heatmap.pdf",
  plot = p_heatmap_top20,
  width = 14,
  height = 18
)
# =========================
# 17. DotPlot for marker genes
# =========================

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
  filename = "cluster_top20_marker_dotplot.pdf",
  plot = p_dot,
  width = 16,
  height = 8
)

# DotPlot for selected marker genes.
DefaultAssay(merged_obj_qc) <- "RNA"

genes_to_plot <- c(
  "SFTPC", "SFTPA1", "SFTPB", "SFTPD", "SLC34A2", "NAPSA", "NKX2-1",
  "AGER", "PDPN", "CLIC5", "CAV1",
  "KRT5", "KRT17",
  "SOX9", "COL1A1", "CDKN2B","EPCAM",
  "MKI67", "TOP2A"
)

p_dot_selected <- DotPlot(
  merged_obj_qc,
  features = genes_to_plot,
  group.by = "seurat_clusters"
) +
  scale_color_gradient(
    low = "grey90",
    high = "purple"
  ) +
  RotatedAxis()

p_dot_selected

ggsave(
  filename = "cluster_marker_selected_dotplot.pdf",
  plot = p_dot_selected,
  width = 16,
  height = 8
)

# DotPlot by experimental condition.
p_dot_selected_by_condition <- DotPlot(
  merged_obj_qc,
  features = genes_to_plot,
  group.by = "condition",
  dot.scale = 8
) +
  RotatedAxis()

p_dot_selected_by_condition

ggsave(
  filename = "marker_selected_dotplot_by_condition.pdf",
  plot = p_dot_selected_by_condition,
  width = 16,
  height = 3
)

# =========================
# 18. Violin plots for gene expression and module scores
# =========================

DefaultAssay(merged_obj_qc) <- "SCT"
Idents(merged_obj_qc) <- "seurat_clusters"

genes_to_check <- c(
  "NKX2-1", "SFTPC", "SFTPA1", "SFTPB",
  "SFTPD", "SLC34A2", "ABCA3", "NAPSA",
  "TM4SF1", "AXIN2", "ETV5",
  "FGFR2", "WIF1", "NKD1",
  "AGER", "PDPN",
  "AQP5", "CAV1", "CLIC5", 
  "KRT5", "KRT17", 
  "MKI67", "TOP2A","UBE2C", "CENPF",
  "KRT8", "KRT7", "CLDN4","KRT19", "KRT17", "MMP7", "LGALS3",
  "ATF3", "FOS", "JUN","HSPA1A", "DDIT3", "CDKN1A"
)

p_vln <- VlnPlot(
  merged_obj_qc,
  features = genes_to_check,
  group.by = "seurat_clusters",
  pt.size = 0,
  ncol = 5
)


ggsave(
  filename = "cluster_gene_expression_violin.pdf",
  plot = p_vln,
  width = 12,
  height = 20
)

###########AT1 module score violin plot


# ============================================================
#Violin plot

p_AT1_signature_violin <- VlnPlot(
  merged_obj_qc,
  features = "AT1_signature1",
  group.by = "condition",
  pt.size = 0
) +
  labs(
    title = "AT1 Signature",
    x = NULL,
    y = "Module Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5
    ),
    # X-axis text
    axis.text.x = element_text(
      size = 30,
      angle = 30,
      hjust = 1
    ),
    # Y-axis text
    axis.text.y = element_text(
      size = 20
    ),
    axis.title.y = element_text(
      size = 25,
      face = "bold"
    ),
    legend.position = "none",
    
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_blank()
  )


ggsave(
  filename = "AT1_signature_violin_by_condition.png",
  plot = p_AT1_signature_violin,
  width = 6,
  height = 7,
  dpi = 300
)
###########AT2 module score violin plot


# ============================================================
#Violin plot

p_AT2_signature_violin <- VlnPlot(
  merged_obj_qc,
  features = "AT2_signature1",
  group.by = "condition",
  pt.size = 0
) +
  labs(
    title = "AT2 Signature",
    x = NULL,
    y = "Module Score"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5
    ),
    # X-axis text
    axis.text.x = element_text(
      size = 30,
      angle = 30,
      hjust = 1
    ),
    # Y-axis text
    axis.text.y = element_text(
      size = 20
    ),
    axis.title.y = element_text(
      size = 25,
      face = "bold"
    ),
    legend.position = "none",
    
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.5
    ),
    axis.line = element_blank()
  )


ggsave(
  filename = "AT2_signature_violin_by_condition.png",
  plot = p_AT2_signature_violin,
  width = 6,
  height = 7,
  dpi = 300
)
# =========================
# 19. AT1 and AT2 signature heatmaps by condition

# AT1 / AT2 signature gene heatmap by condition
# ============================================================
# ============================================================
# AT1 / AT2 signature heatmaps separately
# ============================================================


# ============================================================
# Use RNA assay

DefaultAssay(merged_obj_qc) <- "RNA"

merged_obj_qc <- NormalizeData(
  merged_obj_qc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)
# ============================================================
#  Condition order


# ============================================================
# Check which genes are present

AT1_genes_use <- AT1_genes[
  AT1_genes %in% rownames(
    merged_obj_qc[["RNA"]]
  )
]

AT2_genes_use <- AT2_genes[
  AT2_genes %in% rownames(
    merged_obj_qc[["RNA"]]
  )
]


cat(
  "\nAT1 genes used:",
  length(AT1_genes_use),
  "/",
  length(AT1_genes),
  "\n"
)

cat(
  "Missing AT1 genes:\n"
)

print(
  setdiff(
    AT1_genes,
    AT1_genes_use
  )
)


cat(
  "\nAT2 genes used:",
  length(AT2_genes_use),
  "/",
  length(AT2_genes),
  "\n"
)

cat(
  "Missing AT2 genes:\n"
)

print(
  setdiff(
    AT2_genes,
    AT2_genes_use
  )
)


# ============================================================
#  Average expression by condition


avg_AT1 <- AverageExpression(
  merged_obj_qc,
  assays = "RNA",
  features = AT1_genes_use,
  group.by = "condition",
  layer = "data"
)$RNA

avg_AT2 <- AverageExpression(
  merged_obj_qc,
  assays = "RNA",
  features = AT2_genes_use,
  group.by = "condition",
  layer = "data"
)$RNA

# ============================================================
# Fix condition order

avg_AT1 <- avg_AT1[
  ,
  condition_order,
  drop = FALSE
]

avg_AT2 <- avg_AT2[
  ,
  condition_order,
  drop = FALSE
]


# ============================================================
# Calculate row-wise Z-score
# Each gene is standardized across 4 conditions


AT1_heatmap <- t(
  scale(
    t(avg_AT1)
  )
)

AT2_heatmap <- t(
  scale(
    t(avg_AT2)
  )
)


# ============================================================
# Remove genes with NA / Inf
# Usually genes with no variation across conditions


AT1_heatmap <- AT1_heatmap[
  apply(
    AT1_heatmap,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]


AT2_heatmap <- AT2_heatmap[
  apply(
    AT2_heatmap,
    1,
    function(x) all(is.finite(x))
  ),
  ,
  drop = FALSE
]


AT1_heatmap_sorted <- AT1_heatmap
AT2_heatmap_sorted <- AT2_heatmap

# ============================================================
# Optional:
# cap extreme Z-scores for clearer color visualization
# ============================================================

# AT1_heatmap_sorted[
#   AT1_heatmap_sorted > 2
# ] <- 2
# 
# AT1_heatmap_sorted[
#   AT1_heatmap_sorted < -2
# ] <- -2
# 
# 
# AT2_heatmap_sorted[
#   AT2_heatmap_sorted > 2
# ] <- 2
# 
# AT2_heatmap_sorted[
#   AT2_heatmap_sorted < -2
# ] <- -2


# ============================================================
# Shared color scale

heatmap_colors <- colorRampPalette(
  c(
    "navy",
    "white",
    "firebrick3"
  )
)(100)

heatmap_breaks <- seq(
  -2,
  2,
  length.out = 101
)


# ============================================================
# AT1 heatmap

p_AT1_heatmap <- pheatmap(
  AT1_heatmap_sorted,
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_method = "complete",
  
  color = heatmap_colors,
  breaks = heatmap_breaks,
  
  border_color = NA,
  
  fontsize_row = 7,
  fontsize_col = 11,
  
  angle_col = 90,
  
  main = "AT1 Signature Genes",
  
  scale = "none",
  
  silent = TRUE
)


# ============================================================
# AT2 heatmap


p_AT2_heatmap <- pheatmap(
  AT2_heatmap_sorted,
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_method = "complete",
  
  color = heatmap_colors,
  breaks = heatmap_breaks,
  
  border_color = NA,
  
  fontsize_row = 7,
  fontsize_col = 11,
  
  angle_col = 90,
  
  main = "AT2 Signature Genes",
  
  scale = "none",
  
  silent = TRUE
)


# ============================================================
# Display plots

grid::grid.newpage()
grid::grid.draw(
  p_AT1_heatmap$gtable
)

grid::grid.newpage()
grid::grid.draw(
  p_AT2_heatmap$gtable
)


# ============================================================
# Save AT1 heatmap


png(
  filename = "Heatmap_AT1_signature_clustered.png",
  width = 900,
  height = 3600,
  res = 300
)

grid::grid.draw(
  p_AT1_heatmap$gtable
)
# ============================================================
# Save AT2 heatmap


png(
  filename = "Heatmap_AT2_signature_clustered.png",
  width = 900,
  height = 3600,
  res = 300
)

grid::grid.draw(
  p_AT2_heatmap$gtable
)
# ============================================================
# Save heatmap Z-score data to Excel
# ============================================================
AT1_heatmap_sorted[AT1_heatmap_sorted > 2] <- 2
AT1_heatmap_sorted[AT1_heatmap_sorted < -2] <- -2

AT2_heatmap_sorted[AT2_heatmap_sorted > 2] <- 2
AT2_heatmap_sorted[AT2_heatmap_sorted < -2] <- -2

wb <- createWorkbook()


# ============================================================
# AT1 Z-score sheet
# ============================================================

AT1_zscore_export <- data.frame(
  Gene = rownames(AT1_heatmap_sorted),
  AT1_heatmap_sorted,
  check.names = FALSE
)

addWorksheet(
  wb,
  "AT1_Zscore"
)

writeData(
  wb,
  sheet = "AT1_Zscore",
  x = AT1_zscore_export
)


# ============================================================
# AT2 Z-score sheet
# ============================================================

AT2_zscore_export <- data.frame(
  Gene = rownames(AT2_heatmap_sorted),
  AT2_heatmap_sorted,
  check.names = FALSE
)

addWorksheet(
  wb,
  "AT2_Zscore"
)

writeData(
  wb,
  sheet = "AT2_Zscore",
  x = AT2_zscore_export
)


# ============================================================
# Save Excel file
# ============================================================

saveWorkbook(
  wb,
  file = "AT1_AT2_signature_heatmap_Zscores.xlsx",
  overwrite = TRUE
)
# ============================================================
# 20. Inhibitor downstream gene expression DotPlot
# DMSO / LATSi / MAPKi / LATSi+MAPKi
# ============================================================

DefaultAssay(merged_obj_qc) <- "RNA"


# ============================================================
# Downstream genes
# ============================================================

LATSi_genes <- c(
  "CCN2", #CTGF
  "CCN1", #CYR61
  "ANKRD1",
  "INHBA",
  "TEAD1",
  "TEAD2",
  "TEAD3",
  "TEAD4",
  "WWTR1", #TAZ
  "YAP1"
)

MAPKi_genes <- c(
  "DUSP6",
  "DUSP4",
  "ETV5",
  "ETV4",
  "SPRY4",
  "SPRY2",
  "EGR1",
  "IER3",
  "MAPK1",
  "MAPK3",
  "MAP2K2",
  "MAP2K1"
)

genes_to_plot <- c(
  LATSi_genes,
  MAPKi_genes
)


# ============================================================
# Check genes present in dataset
# ============================================================

genes_present <- genes_to_plot[
  genes_to_plot %in% rownames(merged_obj_qc[["RNA"]])
]

genes_missing <- setdiff(
  genes_to_plot,
  rownames(merged_obj_qc[["RNA"]])
)

cat("\nGenes missing from dataset:\n")
print(genes_missing)


# ============================================================
# Set condition order
# ============================================================


# ============================================================
# DotPlot

# ============================================================
# 21. Combined DotPlot with GSE135893 Control lung reference
# ============================================================
# Compare the four differentiation conditions in the current dataset with
# selected epithelial cell types from the GSE135893 Control lung dataset.
#
# Dot size: percentage of cells expressing each gene
# Dot color: scaled average expression within each dataset
#
# Update this path to the directory containing the GSE135893 files.
gse135893_dir <- "path/to/GSE135893"

# ============================================================
# 21.1 Genes and reference cell types
# ============================================================

combined_genes <- c(
  # AT2
  "SFTPC", "SFTPA1", "SFTPB", "ABCA3", "SLC34A2", "NAPSA",
  
  # Transitional
  "KRT8", "KRT19", "CLDN4", "KRT17", "MMP7", "LGALS3",
  
  # AT1
  "AGER", "PDPN", "CAV1", "CLIC5", "AQP5",
  
  # Basal / airway
  "KRT5", "SCGB3A2", "SCGB1A1",
  
  # Proliferation
  "MKI67", "TOP2A"
)

reference_celltypes <- c(
  "AT2",
  "Transitional AT2",
  "KRT5-/KRT17+",
  "AT1",
  "Basal",
  "SCGB3A2+"
)

# ============================================================
# 21.2 Load GSE135893 reference data
# ============================================================

genes_ref <- read.delim(
  gzfile(file.path(gse135893_dir, "GSE135893_genes.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)
gene_symbols_ref <- genes_ref[[1]]

barcodes_ref <- read.delim(
  gzfile(file.path(gse135893_dir, "GSE135893_barcodes.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)
cell_barcodes_ref <- barcodes_ref[[1]]

meta_ref <- read.csv(
  gzfile(file.path(gse135893_dir, "GSE135893_IPF_metadata.csv.gz")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rownames(meta_ref) <- as.character(meta_ref[[1]])

counts_ref <- Matrix::readMM(
  gzfile(file.path(gse135893_dir, "GSE135893_matrix.mtx.gz"))
)
rownames(counts_ref) <- gene_symbols_ref
colnames(counts_ref) <- cell_barcodes_ref

# ============================================================
# 21.3 Select Control reference cells
# ============================================================

control_cells_ref <- rownames(meta_ref)[
  meta_ref$Diagnosis == "Control" &
    meta_ref$celltype %in% reference_celltypes
]

control_cells_ref <- intersect(
  control_cells_ref,
  colnames(counts_ref)
)

if (length(control_cells_ref) == 0) {
  stop("No requested GSE135893 Control reference cells were found.")
}

# ============================================================
# 21.4 Identify genes shared by both datasets
# ============================================================

DefaultAssay(merged_obj_qc) <- "RNA"

genes_present_ref <- combined_genes[
  combined_genes %in% rownames(counts_ref)
]

genes_present_my <- combined_genes[
  combined_genes %in% rownames(merged_obj_qc[["RNA"]])
]

genes_common <- combined_genes[
  combined_genes %in% genes_present_ref &
    combined_genes %in% genes_present_my
]

genes_missing_ref <- setdiff(
  combined_genes,
  genes_present_ref
)

genes_missing_my <- setdiff(
  combined_genes,
  genes_present_my
)

cat("\nGenes missing from GSE135893:\n")
print(genes_missing_ref)

cat("\nGenes missing from current dataset:\n")
print(genes_missing_my)

cat(
  "\nGenes used in combined DotPlot:",
  length(genes_common),
  "/",
  length(combined_genes),
  "\n"
)

if (length(genes_common) == 0) {
  stop("None of the requested genes are shared by both datasets.")
}

# ============================================================
# 21.5 Extract selected GSE135893 counts
# ============================================================

gene_idx_ref <- match(
  genes_common,
  rownames(counts_ref)
)
gene_idx_ref <- gene_idx_ref[
  !is.na(gene_idx_ref)
]

cell_idx_ref <- match(
  control_cells_ref,
  colnames(counts_ref)
)
cell_idx_ref <- cell_idx_ref[
  !is.na(cell_idx_ref)
]

gene_map_ref <- integer(
  nrow(counts_ref)
)
gene_map_ref[gene_idx_ref] <- seq_along(
  gene_idx_ref
)

cell_map_ref <- integer(
  ncol(counts_ref)
)
cell_map_ref[cell_idx_ref] <- seq_along(
  cell_idx_ref
)

keep_ref <- (
  gene_map_ref[counts_ref@i + 1L] > 0L &
    cell_map_ref[counts_ref@j + 1L] > 0L
)

i_ref <- gene_map_ref[
  counts_ref@i[keep_ref] + 1L
]
j_ref <- cell_map_ref[
  counts_ref@j[keep_ref] + 1L
]
x_ref <- counts_ref@x[
  keep_ref
]

counts_ref_small <- Matrix::sparseMatrix(
  i = i_ref,
  j = j_ref,
  x = x_ref,
  dims = c(
    length(gene_idx_ref),
    length(cell_idx_ref)
  ),
  dimnames = list(
    rownames(counts_ref)[gene_idx_ref],
    colnames(counts_ref)[cell_idx_ref]
  )
)

meta_ref_control <- meta_ref[
  colnames(counts_ref_small),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(meta_ref_control),
    colnames(counts_ref_small)
  )
)

rm(
  counts_ref,
  keep_ref,
  gene_map_ref,
  cell_map_ref,
  i_ref,
  j_ref,
  x_ref
)
gc()

# ============================================================
# 21.6 Normalize GSE135893 reference counts
# ============================================================
# log1p(raw count / total UMI per cell * 10000)

norm_factor_ref <- (
  10000 /
    meta_ref_control$nCount_RNA
)

ref_norm <- counts_ref_small %*%
  Matrix::Diagonal(
    x = norm_factor_ref
  )

ref_norm@x <- log1p(
  ref_norm@x
)

ref_df <- as.data.frame(
  t(
    as.matrix(ref_norm)
  )
)

ref_df$celltype <- as.character(
  meta_ref_control$celltype
)

ref_long <- ref_df %>%
  tidyr::pivot_longer(
    cols = all_of(genes_common),
    names_to = "gene",
    values_to = "expression"
  )

# ============================================================
# 21.7 Calculate reference DotPlot statistics
# ============================================================

ref_dot <- ref_long %>%
  group_by(
    celltype,
    gene
  ) %>%
  summarise(
    avg_expression = mean(
      expression,
      na.rm = TRUE
    ),
    pct_expression = mean(
      expression > 0,
      na.rm = TRUE
    ) * 100,
    .groups = "drop"
  ) %>%
  group_by(gene) %>%
  mutate(
    avg_scaled = as.numeric(
      scale(avg_expression)
    )
  ) %>%
  ungroup()

# Replace NA values produced when a gene has zero variance across groups.
ref_dot$avg_scaled[
  is.na(ref_dot$avg_scaled)
] <- 0

ref_dot$avg_scaled <- pmax(
  pmin(
    ref_dot$avg_scaled,
    2.5
  ),
  -2.5
)

# ============================================================
# 21.8 Calculate DotPlot statistics for the current dataset
# ============================================================

DefaultAssay(merged_obj_qc) <- "RNA"

my_dot <- DotPlot(
  merged_obj_qc,
  features = genes_common,
  group.by = "condition",
  assay = "RNA"
)$data

my_dot2 <- my_dot %>%
  transmute(
    gene = as.character(features.plot),
    group = as.character(id),
    pct_expression = pct.exp,
    avg_scaled = avg.exp.scaled,
    dataset = "Current dataset"
  )

ref_dot2 <- ref_dot %>%
  transmute(
    gene = as.character(gene),
    group = as.character(celltype),
    pct_expression = pct_expression,
    avg_scaled = avg_scaled,
    dataset = "GSE135893 Control"
  )

# ============================================================
# 21.9 Combine current and reference DotPlot data
# ============================================================

combined_dot <- bind_rows(
  my_dot2,
  ref_dot2
)

combined_dot$gene <- factor(
  combined_dot$gene,
  levels = genes_common
)

desired_top_to_bottom <- c(
  "SCGB3A2+",
  "Basal",
  "KRT5-/KRT17+",
  "Transitional AT2",
  "AT2",
  "AT1",
  "LATSi+MAPKi",
  "MAPKi",
  "LATSi",
  "DMSO"
)

# Use numeric y positions to control vertical spacing.
y_spacing <- 0.2

group_positions <- setNames(
  rev(
    seq_along(
      desired_top_to_bottom
    )
  ) * y_spacing,
  desired_top_to_bottom
)

combined_dot$y_pos <- group_positions[
  as.character(
    combined_dot$group
  )
]

if (any(is.na(combined_dot$y_pos))) {
  missing_groups <- unique(
    as.character(
      combined_dot$group
    )[
      is.na(
        combined_dot$y_pos
      )
    ]
  )
  
  stop(
    paste(
      "Unexpected groups in combined DotPlot:",
      paste(
        missing_groups,
        collapse = ", "
      )
    )
  )
}

# ============================================================
# 21.10 Generate combined DotPlot
# ============================================================

p_dot_combined <- ggplot(
  combined_dot,
  aes(
    x = gene,
    y = y_pos
  )
) +
  geom_point(
    aes(
      size = pct_expression,
      color = avg_scaled
    )
  ) +
  scale_size(
    range = c(0, 8),
    limits = c(0, 100)
  ) +
  scale_color_gradientn(
    colors = c(
      "#E5E5E5",
      "#C6DBEF",
      "#6BAED6",
      "#2171B5",
      "#08306B"
    ),
    limits = c(-1, 1.5),
    oob = scales::squish
  ) +
  scale_x_discrete(
    expand = expansion(
      mult = c(0.02, 0.02)
    )
  ) +
  scale_y_continuous(
    breaks = group_positions,
    labels = names(group_positions),
    expand = expansion(
      mult = c(0.03, 0.03)
    )
  ) +
  labs(
    x = NULL,
    y = NULL,
    size = "% expressed",
    color = "Scaled\nexpression"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 20
    ),
    axis.text.y = element_text(
      size = 20
    ),
    axis.title = element_blank(),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.4
    ),
    axis.line = element_blank()
  )

# ============================================================
# 21.11 Save combined DotPlot and numerical data
# ============================================================

ggsave(
  filename = "DotPlot_AT2_to_AT1_current_GSE135893.png",
  plot = p_dot_combined,
  width = 15,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "DotPlot_AT2_to_AT1_current_GSE135893.pdf",
  plot = p_dot_combined,
  width = 10,
  height = 7
)

write.csv(
  combined_dot,
  file = "DotPlot_AT2_to_AT1_current_GSE135893_data.csv",
  row.names = FALSE
)

# ============================================================
# 22. Save final processed object
# ============================================================
# ============================================================

p_inhibitor_dot <- DotPlot(
  merged_obj_qc,
  features = genes_present,
  group.by = "condition",
  dot.scale = 7
) +
  scale_color_gradientn(
    colors = c(    "#2166AC",
                   "#FFFFBF",
                   "#B2182B")
  ) +
  coord_flip() +
  
  labs(
    x = NULL,
    y = NULL,
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "italic"
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    )
  )


ggsave(
  filename = "dot_plot_signaling.png",
  plot = p_inhibitor_dot,
  width = 6,
  height = 8,
  dpi = 300
)
############################ Save final object.
###########################

saveRDS(
  merged_obj_qc,
  file = "SB021_SB024_processed.rds",
  compress = FALSE
)


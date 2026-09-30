# ============================================================
# GSE136831 Control lung
# Growth factor expression across endothelial cell subtypes
# and fibroblasts
#
# Analysis:
#   1. Select Control cells
#   2. Log-normalize raw counts
#   3. Calculate donor-level mean expression
#   4. Calculate cell-type mean expression across donors
#   5. Row-scale expression (z-score) across cell types
#   6. Generate a Control-only heatmap
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(pheatmap)


# ============================================================
# 1. File paths
# ============================================================

# Update this path to the directory containing the GSE136831 files.
gse136831_dir <- "path/to/GSE136831"

count_file <- file.path(
  gse136831_dir,
  "GSE136831_RawCounts_Sparse.mtx.gz"
)

gene_file <- file.path(
  gse136831_dir,
  "GSE136831_AllCells.GeneIDs.txt.gz"
)

barcode_file <- file.path(
  gse136831_dir,
  "GSE136831_AllCells.cellBarcodes.txt.gz"
)

meta_file <- file.path(
  gse136831_dir,
  "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
)


# ============================================================
# 2. Load metadata
# ============================================================

meta <- fread(
  meta_file,
  data.table = FALSE
)


# ============================================================
# 3. Select Control populations
# ============================================================

celltype_order <- c(
  "VE_Capillary_A",
  "VE_Capillary_B",
  "VE_Arterial",
  "VE_Venous",
  "VE_Peribronchial",
  "Fibroblast"
)

meta_sub <- meta %>%
  filter(
    Disease_Identity == "Control",
    Manuscript_Identity %in% celltype_order
  )

if (nrow(meta_sub) == 0) {
  stop("No Control cells were found for the requested cell populations.")
}


# ============================================================
# 4. Growth factor-related genes
# ============================================================

genes_target <- c(
  "AREG",
  "FGF2",
  "EGF",
  "EGFR",
  "FGF4",
  "FGF6",
  "FGF7",
  "FGF10",
  "CSF3",
  "GDNF",
  "CSF2",
  "HBEGF",
  "HGF",
  "IGFBP1",
  "IGFBP2",
  "IGFBP3",
  "IGFBP4",
  "IGFBP6",
  "IGF1",
  "IGF1R",
  "IGF2",
  "CSF1",
  "CSF1R",
  "NGF",
  "NTF3",
  "NTF4",
  "PDGFRA",
  "PDGFRB",
  "PDGFA",
  "PDGFB",
  "PGF",
  "KITLG",
  "KIT",
  "TGFA",
  "TGFB1",
  "TGFB2",
  "TGFB3",
  "VEGFA",
  "KDR",
  "FLT4",
  "FIGF"
)

gene_labels <- c(
  AREG   = "Amphiregulin",
  FGF2   = "bFGF",
  EGF    = "EGF",
  EGFR   = "EGFR",
  FGF4   = "FGF-4",
  FGF6   = "FGF-6",
  FGF7   = "FGF-7",
  FGF10  = "FGF-10",
  CSF3   = "G-CSF",
  GDNF   = "GDNF",
  CSF2   = "GM-CSF",
  HBEGF  = "HB-EGF",
  HGF    = "HGF",
  IGFBP1 = "IGFBP-1",
  IGFBP2 = "IGFBP-2",
  IGFBP3 = "IGFBP-3",
  IGFBP4 = "IGFBP-4",
  IGFBP6 = "IGFBP-6",
  IGF1   = "IGF-I",
  IGF1R  = "IGF-I R",
  IGF2   = "IGF-II",
  CSF1   = "M-CSF",
  CSF1R  = "M-CSF R",
  NGF    = "beta-NGF",
  NTF3   = "NT-3",
  NTF4   = "NT-4",
  PDGFRA = "PDGF Ralpha",
  PDGFRB = "PDGF Rbeta",
  PDGFA  = "PDGF-A",
  PDGFB  = "PDGF-B",
  PGF    = "PLGF",
  KITLG  = "SCF",
  KIT    = "SCF R",
  TGFA   = "TGF-alpha",
  TGFB1  = "TGF-beta1",
  TGFB2  = "TGF-beta2",
  TGFB3  = "TGF-beta3",
  VEGFA  = "VEGF-A",
  KDR    = "VEGF R2",
  FLT4   = "VEGF R3",
  FIGF   = "VEGF-D"
)


# ============================================================
# 5. Load gene IDs and cell barcodes
# ============================================================

genes <- fread(
  gene_file,
  data.table = FALSE
)

barcodes <- fread(
  barcode_file,
  header = FALSE,
  data.table = FALSE
)

barcode_vec <- barcodes[[1]]
gene_names <- genes$HGNC_EnsemblAlt_GeneID


# ============================================================
# 6. Load sparse count matrix
# ============================================================

counts <- Matrix::readMM(
  count_file
)

counts <- as(
  counts,
  "dgCMatrix"
)

if (nrow(counts) != length(gene_names)) {
  stop("The number of matrix rows does not match the number of gene IDs.")
}

if (ncol(counts) != length(barcode_vec)) {
  stop("The number of matrix columns does not match the number of barcodes.")
}

rownames(counts) <- gene_names
colnames(counts) <- barcode_vec


# ============================================================
# 7. Match metadata and count matrix
# ============================================================

cells_keep <- intersect(
  meta_sub$CellBarcode_Identity,
  colnames(counts)
)

if (length(cells_keep) == 0) {
  stop("None of the selected Control cells were found in the count matrix.")
}

meta_sub <- meta_sub %>%
  filter(
    CellBarcode_Identity %in% cells_keep
  )

total_umi <- Matrix::colSums(
  counts
)

total_umi_sub <- total_umi[
  cells_keep
]

if (any(total_umi_sub <= 0)) {
  stop("At least one selected cell has zero total UMI counts.")
}


# ============================================================
# 8. Check target genes
# ============================================================

genes_found <- genes_target[
  genes_target %in% rownames(counts)
]

genes_missing <- setdiff(
  genes_target,
  rownames(counts)
)

cat(
  "\nGrowth factor-related genes found:",
  length(genes_found),
  "/",
  length(genes_target),
  "\n"
)

if (length(genes_missing) > 0) {
  cat("\nGenes missing from GSE136831:\n")
  print(genes_missing)
}

if (length(genes_found) == 0) {
  stop("None of the requested growth factor-related genes were found.")
}

counts_target <- counts[
  genes_found,
  cells_keep,
  drop = FALSE
]


# ============================================================
# 9. Log-normalize expression
# ============================================================
# log1p(raw count / total UMI per cell * 10000)

norm_target <- counts_target

norm_target <- t(
  t(norm_target) /
    total_umi_sub *
    10000
)

norm_target@x <- log1p(
  norm_target@x
)


# ============================================================
# 10. Create expression dataframe
# ============================================================

expr_df <- as.data.frame(
  as.matrix(norm_target)
)

expr_df$Gene <- rownames(
  expr_df
)

expr_long <- expr_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "CellBarcode_Identity",
    values_to = "Expression"
  ) %>%
  left_join(
    meta_sub %>%
      select(
        CellBarcode_Identity,
        Subject_Identity,
        Manuscript_Identity
      ),
    by = "CellBarcode_Identity"
  )

if (any(is.na(expr_long$Subject_Identity))) {
  stop("Metadata could not be matched to all selected cells.")
}


# ============================================================
# 11. Calculate donor-level mean expression
# ============================================================
# Each donor contributes one mean expression value per
# cell type and gene, preventing donors with more cells from
# disproportionately weighting the cell-type average.

donor_expr <- expr_long %>%
  group_by(
    Subject_Identity,
    Manuscript_Identity,
    Gene
  ) %>%
  summarise(
    Expression = mean(
      Expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

donor_expr$Manuscript_Identity <- factor(
  donor_expr$Manuscript_Identity,
  levels = celltype_order
)


# ============================================================
# 12. Calculate cell-type mean expression across donors
# ============================================================

heatmap_df <- donor_expr %>%
  group_by(
    Gene,
    Manuscript_Identity
  ) %>%
  summarise(
    MeanExpression = mean(
      Expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

heatmap_wide <- heatmap_df %>%
  pivot_wider(
    names_from = Manuscript_Identity,
    values_from = MeanExpression
  )

missing_celltypes <- setdiff(
  celltype_order,
  colnames(heatmap_wide)
)

if (length(missing_celltypes) > 0) {
  stop(
    paste(
      "Missing requested cell types:",
      paste(
        missing_celltypes,
        collapse = ", "
      )
    )
  )
}

heatmap_mat <- as.matrix(
  heatmap_wide[
    ,
    celltype_order
  ]
)

rownames(heatmap_mat) <- heatmap_wide$Gene


# ============================================================
# 13. Row-scale expression across cell types
# ============================================================
# Z-scores are calculated independently for each gene across
# the six cell populations. Therefore, the heatmap represents
# relative cell-type enrichment for each gene rather than
# absolute expression differences between different genes.

heatmap_z <- t(
  scale(
    t(heatmap_mat)
  )
)

finite_rows <- apply(
  heatmap_z,
  1,
  function(x) {
    all(
      is.finite(x)
    )
  }
)

genes_removed_zero_variance <- rownames(
  heatmap_z
)[
  !finite_rows
]

if (length(genes_removed_zero_variance) > 0) {
  cat(
    "\nGenes removed because row-wise z-scores could not be calculated:\n"
  )
  print(
    genes_removed_zero_variance
  )
}

heatmap_z <- heatmap_z[
  finite_rows,
  ,
  drop = FALSE
]

if (nrow(heatmap_z) == 0) {
  stop("No genes remained after row-wise z-score filtering.")
}


# ============================================================
# 14. Apply display labels
# ============================================================

display_labels <- gene_labels[
  rownames(heatmap_z)
]

missing_labels <- is.na(
  display_labels
)

display_labels[
  missing_labels
] <- rownames(
  heatmap_z
)[
  missing_labels
]


# ============================================================
# 15. Generate Control-only heatmap
# ============================================================

pheatmap(
  heatmap_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  labels_row = unname(
    display_labels
  ),
  main = "Growth factor expression in control lung",
  filename = "GSE136831_Control_EC_GrowthFactor_Heatmap.png",
  width = 8,
  height = 11
)


# ============================================================
# 16. Export numerical heatmap data
# ============================================================

write.csv(
  heatmap_mat,
  file = "GSE136831_Control_EC_GrowthFactor_mean_expression.csv",
  row.names = TRUE
)

write.csv(
  heatmap_z,
  file = "GSE136831_Control_EC_GrowthFactor_row_zscore.csv",
  row.names = TRUE
)

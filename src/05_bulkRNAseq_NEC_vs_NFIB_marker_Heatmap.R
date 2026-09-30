# ============================================================
# Bulk RNA-seq: NEC vs NFIB marker analysis
# ============================================================
# Dataset: Set1 pooled bulk RNA-seq
# Comparison: NEC vs NFIB
# Excluded sample: 5KDJK6_8
#
# Differential expression:
#   edgeR quasi-likelihood framework with TMM normalization
#   and robust dispersion fitting.
#
# Interpretation:
#   logFC > 0: higher expression in NEC
#   logFC < 0: higher expression in NFIB
#
# Marker expression heatmaps use logCPM calculated before
# DEG filtering so that lower-expression marker genes can still
# be inspected. Marker genes removed by filterByExpr retain
# expression values but have NA DEG statistics.
# ============================================================

# ============================================================
# 0. Packages
# ============================================================

library(edgeR)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)
library(openxlsx)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ============================================================
# 1. File paths
# ============================================================

# Update this path to the directory containing the bulk RNA-seq input files.
bulk_data_dir <- "path/to/bulkRNAseq"

# ============================================================
# 2. User settings
# ============================================================

rdata_file <- file.path(
  bulk_data_dir,
  "Set1_pooled_edgeR_analysis.RData"
)

metadata_file <- file.path(
  bulk_data_dir,
  "metadata_set1(1)(1).csv"
)

exclude_samples <- c(
  "5KDJK6_8"
)

group_order <- c(
  "NEC",
  "NFIB"
)


# ============================================================
# 3. Load previous bulk RNA-seq data
# ============================================================

load(
  rdata_file
)

# ============================================================
# 4. Check dge_set1
# ============================================================
#
# We need raw counts to truly recalculate DEG after
# removing 5KDJK6_8.
# ============================================================

if (!exists("dge_set1")) {
  
  stop(
    paste0(
      "\n'dge_set1' was not found in ",
      rdata_file,
      ".\n",
      "To remove 5KDJK6_8 from the DEG statistics, ",
      "the raw-count DGEList object dge_set1 is required."
    )
  )
}


if (is.null(dge_set1$counts)) {
  
  stop(
    "dge_set1 exists, but dge_set1$counts was not found."
  )
}


# ============================================================
# 5. Read metadata
# ============================================================

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_cols <- c(
  "sample",
  "condition"
)


missing_cols <- setdiff(
  required_cols,
  colnames(metadata)
)


if (length(missing_cols) > 0) {
  
  stop(
    paste0(
      "Missing metadata columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 6. Keep NEC / NFIB and remove 5KDJK6_8
# ============================================================

meta_plot <- metadata %>%
  
  dplyr::filter(
    condition %in% c(
      "NEC",
      "NFIB"
    )
  ) %>%
  
  dplyr::filter(
    !sample %in% exclude_samples
  ) %>%
  
  dplyr::mutate(
    
    condition = factor(
      condition,
      levels = c(
        "NEC",
        "NFIB"
      )
    )
    
  ) %>%
  
  dplyr::arrange(
    condition,
    sample
  )


# Confirm excluded sample is gone

if (
  any(
    meta_plot$sample %in% exclude_samples
  )
) {
  
  stop(
    "Excluded sample is still present in metadata."
  )
}


# ============================================================
# 7. Check count-matrix sample names
# ============================================================

missing_samples <- setdiff(
  meta_plot$sample,
  colnames(
    dge_set1$counts
  )
)


if (length(missing_samples) > 0) {
  
  
  
  stop(
    "Sample matching failed."
  )
}


# ============================================================
# 8. Create NEC / NFIB raw-count matrix
# ============================================================

counts_2group <- dge_set1$counts[
  ,
  meta_plot$sample,
  drop = FALSE
]


# ============================================================
# 9. Create DGEList
# ============================================================
#
# Important:
#
# For DEG analysis, set reference = NFIB.
# Therefore:
#
# logFC > 0 = higher in NEC
# logFC < 0 = higher in NFIB
#
# This preserves the interpretation of your previous code.
# ============================================================

dge_all <- edgeR::DGEList(
  
  counts = counts_2group,
  
  group = factor(
    meta_plot$condition,
    levels = c(
      "NFIB",
      "NEC"
    )
  )
  
)


# ============================================================
# 10. TMM normalization
# ============================================================

dge_all <- edgeR::calcNormFactors(
  dge_all
)


# ============================================================
# 11. Calculate logCPM for ALL genes
# ============================================================
#
# Used for marker-expression heatmaps.
# This is calculated before DEG filtering so that
# lower-expression markers can still be inspected.
# ============================================================

logCPM <- edgeR::cpm(
  dge_all,
  log = TRUE,
  prior.count = 2
)


# ============================================================
# 12. Filter genes for DEG analysis
# ============================================================

keep_deg <- edgeR::filterByExpr(
  dge_all
)


dge_deg <- dge_all[
  keep_deg,
  ,
  keep.lib.sizes = FALSE
]


dge_deg <- edgeR::calcNormFactors(
  dge_deg
)


# ============================================================
# 13. Create design matrix
# ============================================================

group_deg <- factor(
  meta_plot$condition,
  levels = c(
    "NFIB",
    "NEC"
  )
)


design <- model.matrix(
  ~ group_deg
)


rownames(
  design
) <- meta_plot$sample


# ============================================================
# 14. Estimate dispersion
# ============================================================

dge_deg <- edgeR::estimateDisp(
  dge_deg,
  design
)


# ============================================================
# 15. Fit edgeR QL model
# ============================================================

fit_NEC_NFIB <- edgeR::glmQLFit(
  dge_deg,
  design,
  robust = TRUE
)


# ============================================================
# 16. NEC vs NFIB test
# ============================================================
#
# coefficient group_degNEC:
#
# positive logFC = NEC > NFIB
# negative logFC = NFIB > NEC
# ============================================================

qlf_NEC_NFIB <- edgeR::glmQLFTest(
  fit_NEC_NFIB,
  coef = "group_degNEC"
)


deg_NEC_NFIB <- edgeR::topTags(
  
  qlf_NEC_NFIB,
  
  n = Inf,
  
  sort.by = "none"
  
)$table


deg_NEC_NFIB$gene_id <- rownames(
  deg_NEC_NFIB
)


deg_NEC_NFIB <- deg_NEC_NFIB %>%
  
  dplyr::select(
    gene_id,
    everything()
  )


# ============================================================
# 17. Save recalculated full DEG table
# ============================================================

openxlsx::write.xlsx(
  
  deg_NEC_NFIB,
  
  file =
    "NEC_vs_NFIB_edgeR_DEG_excluding_5KDJK6_8.xlsx",
  
  overwrite = TRUE
  
)


# ============================================================
# 18. Define ALL marker genes
# ============================================================

marker_list <- list(
  
  Canonical_Endothelial = c(
    "PECAM1",
    "VWF",
    "CDH5",
    "KDR"
  ),
  
  
  Canonical_Fibroblast = c(
    "COL1A1",
    "COL1A2",
    "PDGFRA",
    "DCN",
    "VIM"
  ),
  
  
  VE_Capillary_A = c(
    "HPGD",
    "CA4",
    "SOSTDC1",
    "IL1RL1",
    "EDNRB"
  ),
  
  
  VE_Capillary_B = c(
    "FCN3",
    "IL7R",
    "BTNL9",
    "SLC6A4",
    "SEMA3G"
  ),
  
  
  VE_Arterial = c(
    "DKK2",
    "HEY1",
    "GJA5",
    "IGFBP3",
    "ARL15"
  ),
  
  
  VE_Venous = c(
    "ACKR1",
    "SULT1E1",
    "HDAC9",
    "ADGRG6",
    "FAM155A"
  ),
  
  
  VE_Peribronchial = c(
    "PLVAP",
    "SELE",
    "COL15A1",
    "MPZL2",
    "POSTN"
  )
)


# ============================================================
# 19. Define marker-group order
# ============================================================

subtype_order <- c(
  "Canonical_Endothelial",
  "Canonical_Fibroblast",
  "VE_Capillary_A",
  "VE_Capillary_B",
  "VE_Arterial",
  "VE_Venous",
  "VE_Peribronchial"
)


gene_order <- unlist(
  marker_list,
  use.names = FALSE
)


# ============================================================
# 20. Marker table
# ============================================================

marker_df <- dplyr::bind_rows(
  
  lapply(
    
    names(
      marker_list
    ),
    
    function(subtype) {
      
      data.frame(
        
        Subtype =
          subtype,
        
        Gene =
          marker_list[[subtype]],
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
  )
)


marker_df <- marker_df %>%
  
  dplyr::mutate(
    
    Subtype = factor(
      Subtype,
      levels = subtype_order
    ),
    
    Gene = factor(
      Gene,
      levels = gene_order
    )
    
  ) %>%
  
  dplyr::arrange(
    Subtype,
    Gene
  )


# ============================================================
# 21. Convert gene symbols to Ensembl IDs
# ============================================================

gene_map <- AnnotationDbi::select(
  
  org.Hs.eg.db,
  
  keys =
    unique(
      as.character(
        marker_df$Gene
      )
    ),
  
  keytype = "SYMBOL",
  
  columns = c(
    "SYMBOL",
    "ENSEMBL"
  )
  
)


gene_map <- gene_map %>%
  
  dplyr::filter(
    !is.na(
      ENSEMBL
    )
  ) %>%
  
  dplyr::distinct(
    SYMBOL,
    ENSEMBL
  )


# ============================================================
# 22. Join marker table and mapping
# ============================================================

marker_map <- marker_df %>%
  
  dplyr::mutate(
    Gene = as.character(
      Gene
    )
  ) %>%
  
  dplyr::left_join(
    
    gene_map,
    
    by = c(
      "Gene" = "SYMBOL"
    )
    
  ) %>%
  
  dplyr::mutate(
    
    Found_in_bulk =
      ENSEMBL %in%
      rownames(
        logCPM
      ),
    
    Found_in_DEG =
      ENSEMBL %in%
      deg_NEC_NFIB$gene_id
    
  )


# ============================================================
# 23. Handle multiple Ensembl IDs
# ============================================================
#
# If one gene symbol maps to >1 Ensembl IDs,
# retain the Ensembl ID with the highest mean logCPM.
# ============================================================

marker_map_use <- marker_map %>%
  
  dplyr::filter(
    !is.na(
      ENSEMBL
    ),
    Found_in_bulk
  ) %>%
  
  dplyr::rowwise() %>%
  
  dplyr::mutate(
    
    Mean_logCPM =
      mean(
        
        logCPM[
          ENSEMBL,
          meta_plot$sample,
          drop = TRUE
        ],
        
        na.rm = TRUE
        
      )
    
  ) %>%
  
  dplyr::ungroup() %>%
  
  dplyr::group_by(
    Subtype,
    Gene
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      Mean_logCPM
    ),
    .by_group = TRUE
  ) %>%
  
  dplyr::slice_head(
    n = 1
  ) %>%
  
  dplyr::ungroup()


marker_map_use <- marker_map_use %>%
  
  dplyr::mutate(
    
    Subtype = factor(
      Subtype,
      levels = subtype_order
    ),
    
    Gene = factor(
      Gene,
      levels = gene_order
    )
    
  ) %>%
  
  dplyr::arrange(
    Subtype,
    Gene
  )


# ============================================================
# 24. Extract marker expression matrix
# ============================================================

expr_mat <- logCPM[
  
  marker_map_use$ENSEMBL,
  
  meta_plot$sample,
  
  drop = FALSE
  
]


rownames(
  expr_mat
) <- as.character(
  marker_map_use$Gene
)


row_subtype <- as.character(
  marker_map_use$Subtype
)


# ============================================================
# 25. Calculate NEC / NFIB group mean logCPM
# ============================================================

group_mean_mat <- sapply(
  
  group_order,
  
  function(grp) {
    
    grp_samples <- meta_plot$sample[
      meta_plot$condition == grp
    ]
    
    
    rowMeans(
      
      expr_mat[
        ,
        grp_samples,
        drop = FALSE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
)


colnames(
  group_mean_mat
) <- group_order


# ============================================================
# 26. Row gaps
# ============================================================

subtype_counts <- table(
  
  factor(
    row_subtype,
    levels = subtype_order
  )
  
)


subtype_counts_nonzero <- subtype_counts[
  subtype_counts > 0
]


gaps_row <- cumsum(
  subtype_counts_nonzero
)


if (
  length(
    gaps_row
  ) > 1
) {
  
  gaps_row <- gaps_row[
    -length(
      gaps_row
    )
  ]
  
} else {
  
  gaps_row <- NULL
  
}


# ============================================================
# 27. Group-mean row Z-score
# ============================================================

gene_sd <- apply(
  group_mean_mat,
  1,
  sd,
  na.rm = TRUE
)


keep_var <- is.finite(
  gene_sd
) &
  gene_sd > 0


group_mean_z_mat <- group_mean_mat[
  keep_var,
  ,
  drop = FALSE
]


row_subtype_z <- row_subtype[
  keep_var
]


group_mean_z_mat <- t(
  
  scale(
    
    t(
      group_mean_z_mat
    )
    
  )
  
)


group_mean_z_plot <- group_mean_z_mat


group_mean_z_plot[
  group_mean_z_plot > 2
] <- 2


group_mean_z_plot[
  group_mean_z_plot < -2
] <- -2


subtype_counts_z <- table(
  
  factor(
    row_subtype_z,
    levels = subtype_order
  )
  
)


subtype_counts_z <- subtype_counts_z[
  subtype_counts_z > 0
]


gaps_row_z <- cumsum(
  subtype_counts_z
)


if (
  length(
    gaps_row_z
  ) > 1
) {
  
  gaps_row_z <- gaps_row_z[
    -length(
      gaps_row_z
    )
  ]
  
} else {
  
  gaps_row_z <- NULL
  
}


# ============================================================
# 28. Output directory
# ============================================================

output_dir <- file.path(
  bulk_data_dir,
  "NEC_vs_NFIB_markers_excluding_5KDJK6_8"
)


dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 29. Group mean logCPM heatmap
# ============================================================

group_logcpm_png <- file.path(
  
  output_dir,
  
  "NEC_vs_NFIB_group_mean_logCPM_heatmap.png"
  
)


pheatmap::pheatmap(
  
  group_mean_mat,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  gaps_row = gaps_row,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 12,
  
  border_color = NA,
  
  angle_col = 0,
  
  main =
    "Endothelial/Fibroblast markers: NEC vs NFIB mean logCPM",
  
  filename =
    group_logcpm_png,
  
  width = 5.5,
  
  height = 10
  
)


# ============================================================
# 30. Group mean Z-score heatmap
# ============================================================

z_breaks <- seq(
  -2,
  2,
  length.out = 101
)


z_colors <- colorRampPalette(
  
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
  
)(
  100
)


group_z_png <- file.path(
  
  output_dir,
  
  "NEC_vs_NFIB_group_mean_Zscore_heatmap.png"
  
)


pheatmap::pheatmap(
  
  group_mean_z_plot,
  
  color = z_colors,
  
  breaks = z_breaks,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  gaps_row = gaps_row_z,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 12,
  
  border_color = NA,
  
  angle_col = 0,
  
  main =
    "Endothelial/Fibroblast markers: NEC vs NFIB mean row Z-score",
  
  filename =
    group_z_png,
  
  width = 5.5,
  
  height = 10
  
)


# ============================================================
# 31. Determine subject column
# ============================================================

possible_subject_cols <- c(
  "subject",
  "Subject",
  "subject_id",
  "Subject_ID",
  "subjectID",
  "SubjectID"
)


subject_col <- possible_subject_cols[
  possible_subject_cols %in%
    colnames(
      meta_plot
    )
]


if (
  length(
    subject_col
  ) > 0
) {
  
  subject_col <- subject_col[
    1
  ]
  
  
  meta_plot$SubjectLabel <- as.character(
    meta_plot[
      [
        subject_col
      ]
    ]
  )
  
  
  
  
} else {
  
  meta_plot$SubjectLabel <- as.character(
    meta_plot$sample
  )
  
  
  
}


# ============================================================
# 32. Subject display labels
# ============================================================

meta_plot <- meta_plot %>%
  
  dplyr::mutate(
    
    DisplaySubject =
      paste0(
        
        as.character(
          condition
        ),
        
        "_",
        
        SubjectLabel
        
      )
    
  ) %>%
  
  dplyr::arrange(
    condition,
    SubjectLabel
  )


# ============================================================
# 33. Subject-level expression matrix
# ============================================================

subject_sample_order <- meta_plot$sample


subject_expr <- expr_mat[
  ,
  subject_sample_order,
  drop = FALSE
]


colnames(
  subject_expr
) <- meta_plot$DisplaySubject


# ============================================================
# 34. Column annotation
# ============================================================

annotation_subject <- data.frame(
  
  Group =
    as.character(
      meta_plot$condition
    ),
  
  Subject =
    meta_plot$SubjectLabel
  
)


rownames(
  annotation_subject
) <- meta_plot$DisplaySubject


annotation_subject <- annotation_subject[
  colnames(
    subject_expr
  ),
  ,
  drop = FALSE
]


stopifnot(
  
  identical(
    
    rownames(
      annotation_subject
    ),
    
    colnames(
      subject_expr
    )
    
  )
  
)


# ============================================================
# 35. Gap between NEC and NFIB
# ============================================================

subject_group_counts <- table(
  
  factor(
    
    meta_plot$condition,
    
    levels = c(
      "NEC",
      "NFIB"
    )
    
  )
  
)


gaps_col_subject <- cumsum(
  subject_group_counts
)


if (
  length(
    gaps_col_subject
  ) > 1
) {
  
  gaps_col_subject <- gaps_col_subject[
    -length(
      gaps_col_subject
    )
  ]
  
} else {
  
  gaps_col_subject <- NULL
  
}


# ============================================================
# 36. Subject-level logCPM heatmap
# ============================================================

subject_logcpm_png <- file.path(
  
  output_dir,
  
  "NEC_vs_NFIB_subject_logCPM_heatmap.png"
  
)


pheatmap::pheatmap(
  
  subject_expr,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  gaps_row = gaps_row,
  
  gaps_col = gaps_col_subject,
  
  annotation_col =
    annotation_subject,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 8,
  
  angle_col = 45,
  
  border_color = NA,
  
  main =
    "Endothelial/Fibroblast markers: subject-level logCPM",
  
  filename =
    subject_logcpm_png,
  
  width = 10,
  
  height = 10
  
)


# ============================================================
# 37. Subject-level Z-score
# ============================================================

subject_gene_sd <- apply(
  subject_expr,
  1,
  sd,
  na.rm = TRUE
)


keep_subject_var <- is.finite(
  subject_gene_sd
) &
  subject_gene_sd > 0


subject_expr_z_use <- subject_expr[
  keep_subject_var,
  ,
  drop = FALSE
]


row_subtype_subject_z <- row_subtype[
  keep_subject_var
]


subject_z <- t(
  
  scale(
    
    t(
      subject_expr_z_use
    )
    
  )
  
)


subject_z_plot <- subject_z


subject_z_plot[
  subject_z_plot > 2
] <- 2


subject_z_plot[
  subject_z_plot < -2
] <- -2


subject_subtype_counts_z <- table(
  
  factor(
    
    row_subtype_subject_z,
    
    levels = subtype_order
    
  )
  
)


subject_subtype_counts_z <- subject_subtype_counts_z[
  subject_subtype_counts_z > 0
]


gaps_row_subject_z <- cumsum(
  subject_subtype_counts_z
)


if (
  length(
    gaps_row_subject_z
  ) > 1
) {
  
  gaps_row_subject_z <- gaps_row_subject_z[
    -length(
      gaps_row_subject_z
    )
  ]
  
} else {
  
  gaps_row_subject_z <- NULL
  
}


# ============================================================
# 38. Subject-level Z-score heatmap
# ============================================================

subject_z_png <- file.path(
  
  output_dir,
  
  "NEC_vs_NFIB_subject_Zscore_heatmap.png"
  
)


pheatmap::pheatmap(
  
  subject_z_plot,
  
  color = z_colors,
  
  breaks = z_breaks,
  
  cluster_rows = FALSE,
  
  cluster_cols = FALSE,
  
  gaps_row =
    gaps_row_subject_z,
  
  gaps_col =
    gaps_col_subject,
  
  annotation_col =
    annotation_subject,
  
  show_rownames = TRUE,
  
  show_colnames = TRUE,
  
  fontsize_row = 9,
  
  fontsize_col = 8,
  
  angle_col = 45,
  
  border_color = NA,
  
  main =
    "Endothelial/Fibroblast markers: subject-level row Z-score",
  
  filename =
    subject_z_png,
  
  width = 10,
  
  height = 10
  
)


# ============================================================
# 39. Build marker statistics table
# ============================================================
#
# Keep marker genes even if they were removed by filterByExpr.
# In that case DEG statistics will be NA.
# ============================================================

marker_stats <- marker_map_use %>%
  
  dplyr::mutate(
    
    Gene =
      as.character(
        Gene
      ),
    
    Subtype =
      as.character(
        Subtype
      )
    
  ) %>%
  
  dplyr::left_join(
    
    deg_NEC_NFIB %>%
      
      dplyr::select(
        gene_id,
        logFC,
        logCPM,
        F,
        PValue,
        FDR
      ),
    
    by = c(
      "ENSEMBL" = "gene_id"
    )
    
  )


# ============================================================
# 40. NEC / NFIB mean logCPM
# ============================================================

nec_samples <- meta_plot$sample[
  meta_plot$condition == "NEC"
]


nfib_samples <- meta_plot$sample[
  meta_plot$condition == "NFIB"
]


# Ensure 5KDJK6_8 is NOT there

if (
  "5KDJK6_8" %in%
  nfib_samples
) {
  
  stop(
    "5KDJK6_8 is unexpectedly still included."
  )
}


marker_stats$NEC_mean_logCPM <- sapply(
  
  marker_stats$ENSEMBL,
  
  function(gene_id) {
    
    mean(
      
      logCPM[
        gene_id,
        nec_samples,
        drop = TRUE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
)


marker_stats$NFIB_mean_logCPM <- sapply(
  
  marker_stats$ENSEMBL,
  
  function(gene_id) {
    
    mean(
      
      logCPM[
        gene_id,
        nfib_samples,
        drop = TRUE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
)


# ============================================================
# 41. Mean-expression difference
# ============================================================

marker_stats <- marker_stats %>%
  
  dplyr::mutate(
    
    Mean_logCPM_difference =
      NEC_mean_logCPM -
      NFIB_mean_logCPM
    
  )


# ============================================================
# 42. Significance / direction
# ============================================================

marker_stats <- marker_stats %>%
  
  dplyr::mutate(
    
    Status =
      dplyr::case_when(
        
        !is.na(
          FDR
        ) &
          FDR < 0.05 &
          logFC > 0
        ~ "NEC_high_FDR<0.05",
        
        
        !is.na(
          FDR
        ) &
          FDR < 0.05 &
          logFC < 0
        ~ "NFIB_high_FDR<0.05",
        
        
        TRUE
        ~ "NS"
        
      )
    
  )


# ============================================================
# 43. Order marker statistics
# ============================================================

marker_stats <- marker_stats %>%
  
  dplyr::mutate(
    
    Subtype =
      factor(
        Subtype,
        levels = subtype_order
      ),
    
    Gene =
      factor(
        Gene,
        levels = gene_order
      )
    
  ) %>%
  
  dplyr::arrange(
    Subtype,
    Gene
  )


# ============================================================
# 44. Print all marker results
# ============================================================


# ============================================================
# 45. Significant NEC-high markers
# ============================================================

NEC_high_markers <- marker_stats %>%
  
  dplyr::filter(
    
    !is.na(
      FDR
    ),
    
    FDR < 0.05,
    
    logFC > 0
    
  ) %>%
  
  dplyr::arrange(
    Subtype,
    FDR
  )


# ============================================================
# 46. Significant NFIB-high markers
# ============================================================

NFIB_high_markers <- marker_stats %>%
  
  dplyr::filter(
    
    !is.na(
      FDR
    ),
    
    FDR < 0.05,
    
    logFC < 0
    
  ) %>%
  
  dplyr::arrange(
    Subtype,
    FDR
  )


# ============================================================
# 47. Marker-category summary
# ============================================================

subtype_summary <- marker_stats %>%
  
  dplyr::group_by(
    Subtype
  ) %>%
  
  dplyr::summarise(
    
    Markers_total =
      dplyr::n(),
    
    
    Markers_in_DEG =
      sum(
        !is.na(
          FDR
        )
      ),
    
    
    NEC_high_significant =
      sum(
        
        !is.na(
          FDR
        ) &
          FDR < 0.05 &
          logFC > 0,
        
        na.rm = TRUE
        
      ),
    
    
    NFIB_high_significant =
      sum(
        
        !is.na(
          FDR
        ) &
          FDR < 0.05 &
          logFC < 0,
        
        na.rm = TRUE
        
      ),
    
    
    Mean_logFC =
      mean(
        logFC,
        na.rm = TRUE
      ),
    
    
    Median_logFC =
      median(
        logFC,
        na.rm = TRUE
      ),
    
    
    .groups = "drop"
    
  )


# ============================================================
# 48. Save marker statistics Excel
# ============================================================

statistics_xlsx <- file.path(
  
  output_dir,
  
  "Endothelial_Fibroblast_markers_NEC_vs_NFIB_statistics.xlsx"
  
)


openxlsx::write.xlsx(
  
  list(
    
    Samples_used =
      as.data.frame(
        meta_plot
      ),
    
    
    All_markers =
      as.data.frame(
        marker_stats
      ),
    
    
    NEC_high_FDR005 =
      as.data.frame(
        NEC_high_markers
      ),
    
    
    NFIB_high_FDR005 =
      as.data.frame(
        NFIB_high_markers
      ),
    
    
    Category_summary =
      as.data.frame(
        subtype_summary
      ),
    
    
    Full_DEG =
      as.data.frame(
        deg_NEC_NFIB
      )
    
  ),
  
  file =
    statistics_xlsx,
  
  overwrite = TRUE
  
)


# ============================================================
# 49. Marker logFC plot
# ============================================================

plot_df <- marker_stats %>%
  
  dplyr::filter(
    !is.na(
      logFC
    )
  ) %>%
  
  dplyr::mutate(
    
    Gene =
      factor(
        
        as.character(
          Gene
        ),
        
        levels =
          rev(
            
            gene_order[
              gene_order %in%
                as.character(
                  Gene
                )
            ]
            
          )
        
      )
    
  )


p_logFC <- ggplot(
  
  plot_df,
  
  aes(
    x = logFC,
    y = Gene
  )
  
) +
  
  geom_point(
    
    aes(
      shape = Status
    ),
    
    size = 3
    
  ) +
  
  geom_vline(
    
    xintercept = 0,
    
    linetype = "dashed"
    
  ) +
  
  facet_grid(
    
    Subtype ~ .,
    
    scales = "free_y",
    
    space = "free_y"
    
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  labs(
    
    title =
      "Endothelial/Fibroblast markers: NEC vs NFIB",
    
    subtitle =
      "5KDJK6_8 excluded",
    
    x =
      "log2 fold change (NEC / NFIB)",
    
    y = NULL
    
  )


logfc_png <- file.path(
  
  output_dir,
  
  "Endothelial_Fibroblast_markers_NEC_vs_NFIB_logFC.png"
  
)


ggsave(
  
  filename =
    logfc_png,
  
  plot =
    p_logFC,
  
  width = 8,
  
  height = 11,
  
  dpi = 300
  
)

# ============================================================
# 50. Save analysis objects
# ============================================================

save(
  
  dge_all,
  dge_deg,
  fit_NEC_NFIB,
  qlf_NEC_NFIB,
  deg_NEC_NFIB,
  logCPM,
  marker_stats,
  meta_plot,
  
  file = file.path(
    
    output_dir,
    
    "NEC_vs_NFIB_analysis_excluding_5KDJK6_8.RData"
    
  )
  
)

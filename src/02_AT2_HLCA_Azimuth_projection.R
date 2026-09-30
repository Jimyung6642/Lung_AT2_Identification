# Azimuth projection of primary and cultured human AT2 cells onto the HLCA Human Lung v2 reference
# Samples: SB001 (P0), SB002 (P5), SB003 (P12)


# ============================================================
# 0. Packages and working directory
# ============================================================
library(Seurat)
library(SeuratData)
library(Azimuth)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(ggrepel)

# Set the working directory to the Azimuth analysis directory before running this script.
# setwd("path/to/Azimuth")
passage_levels <- c("P00", "P05", "P12")

# ============================================================
# 1. Load and clean query object
# ============================================================
merged_obj_qc <- readRDS("../SB001_SB003_processed.rds")

merged_obj_qc <- JoinLayers(
  object = merged_obj_qc,
  assay = "RNA"
)

counts_mat <- GetAssayData(
  object = merged_obj_qc,
  assay = "RNA",
  layer = "counts"
)

meta_df <- merged_obj_qc@meta.data
stopifnot(identical(colnames(counts_mat), rownames(meta_df)))

options(Seurat.object.assay.version = "v3")

query_clean <- CreateSeuratObject(
  counts = counts_mat,
  meta.data = meta_df,
  project = "P0_P5_P12_AT2"
)

DefaultAssay(query_clean) <- "RNA"

# ============================================================
# 2. Passage metadata
# ============================================================
if (!"passage" %in% colnames(query_clean@meta.data)) {
  query_clean$passage <- query_clean$orig.ident
}

query_clean$passage <- as.character(query_clean$passage)
query_clean$passage[query_clean$passage == "SB001"] <- "P00"
query_clean$passage[query_clean$passage == "SB002"] <- "P05"
query_clean$passage[query_clean$passage == "SB003"] <- "P12"

query_clean$passage <- factor(
  query_clean$passage,
  levels = passage_levels
)


# ============================================================
# 3. Run Azimuth
# ============================================================
query_az_hlca <- RunAzimuth(
  query = query_clean,
  reference = "lungref"
)


required_cols <- c(
  "predicted.ann_level_1",
  "predicted.ann_level_3"
)

missing_cols <- setdiff(
  required_cols,
  colnames(query_az_hlca@meta.data)
)

if (length(missing_cols) > 0) {
  stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
}

# ============================================================
# 4. Level-1 epithelial query subset
# ============================================================

query_az_epi <- subset(
  query_az_hlca,
  subset = predicted.ann_level_1 == "Epithelial"
)

query_az_epi$passage <- factor(
  query_az_epi$passage,
  levels = passage_levels
)

# ============================================================
# 5. Load HLCA lung reference and extract epithelial coordinates
# ============================================================
# ============================================================
# Load Azimuth lung reference
# ============================================================
lungref_path <- system.file(package = "lungref.SeuratData")

if (lungref_path == "") {
  stop("Install lungref first with InstallData('lungref').")
}

ref_rds_path <- list.files(
  lungref_path,
  pattern = "^ref\\.Rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(ref_rds_path) == 0) {
  stop("Could not find ref.Rds.")
}

reference_dir <- dirname(ref_rds_path[1])
lung_reference <- LoadReference(path = reference_dir)

# ============================================================
# Extract level-1 epithelial reference coordinates
# ============================================================
ref_df <- as.data.frame(
  Embeddings(lung_reference$plot, reduction = "refUMAP")
)

colnames(ref_df)[1:2] <- c("UMAP_1", "UMAP_2")
ref_df$cell <- rownames(ref_df)

ref_meta <- lung_reference$plot@meta.data
ref_df$ann_level_1 <- ref_meta[ref_df$cell, "ann_level_1"]
ref_df$ann_level_3 <- ref_meta[ref_df$cell, "ann_level_3"]

ref_epi_df <- ref_df %>%
  filter(ann_level_1 == "Epithelial")

# ============================================================
# Extract epithelial query coordinates
# ============================================================
query_epi_df <- as.data.frame(
  Embeddings(query_az_epi, reduction = "ref.umap")
)

colnames(query_epi_df)[1:2] <- c("UMAP_1", "UMAP_2")
query_epi_df$cell <- rownames(query_epi_df)

query_epi_meta <- query_az_epi@meta.data
query_epi_df$passage <- query_epi_meta[query_epi_df$cell, "passage"]
query_epi_df$predicted_ann_level_3 <- query_epi_meta[
  query_epi_df$cell,
  "predicted.ann_level_3"
]

query_epi_df$passage <- factor(
  query_epi_df$passage,
  levels = passage_levels
)

# Define reference-based cell-type order and colors.

if (is.factor(ref_epi_df$ann_level_3)) {
  ref_levels <- levels(droplevels(ref_epi_df$ann_level_3))
} else {
  ref_levels <- unique(ref_epi_df$ann_level_3[!is.na(ref_epi_df$ann_level_3)])
}

ref_levels <- ref_levels[!is.na(ref_levels)]

# Same hue palette style as ggplot default
ref_colors <- setNames(
  scales::hue_pal()(length(ref_levels)),
  ref_levels
)

# ============================================================
# 6. Level-3 epithelial composition by passage
# ============================================================
level3_count <- table(
  Passage = query_az_epi$passage,
  Cell_type = query_az_epi$predicted.ann_level_3
)

level3_percent <- prop.table(
  level3_count,
  margin = 1
) * 100

level3_percent_df <- as.data.frame(level3_percent)

colnames(level3_percent_df) <- c(
  "Passage",
  "Cell_type",
  "Percent"
)

# Add missing reference cell types as 0% to retain a complete legend.


level3_percent_df <- level3_percent_df %>%
  tidyr::complete(
    Passage = passage_levels,
    Cell_type = ref_levels,
    fill = list(
      Percent = 0
    )
  )


# Set plotting order.

level3_percent_df$Passage <- factor(
  level3_percent_df$Passage,
  levels = passage_levels
)

level3_percent_df$Cell_type <- factor(
  level3_percent_df$Cell_type,
  levels = ref_levels
)

# Save level-3 counts and percentages.


write.xlsx(
  list(
    Level3_counts = as.data.frame(level3_count),
    Level3_percent = level3_percent_df
  ),
  file = "Azimuth_epithelial_ann_level3_composition.xlsx",
  overwrite = TRUE
)


# Plot level-3 epithelial composition by passage.


p_level3_bar <- ggplot(
  level3_percent_df,
  aes(
    x = Passage,
    y = Percent,
    fill = Cell_type
  )
) +
  
  geom_col(
    width = 0.75
  ) +
  
  scale_fill_manual(
    values = ref_colors,
    limits = ref_levels,
    breaks = ref_levels,
    drop = FALSE,
    name = "Cell subset"
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  
  scale_x_discrete(
    labels = c(
      P00 = "Passage 0",
      P05 = "Passage 5",
      P12 = "Passage 12"
    )
  ) +  
  theme_classic() +
  
  theme(
    # Axis text
    axis.text.x = element_text(
      size = 20,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 18
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 35,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 35,
      face = "bold"
    ),
    
    # Legend
    legend.title = element_text(
      size = 25,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 25
    )
  ) +
  
  labs(
    x = NULL,
    y = "Percent of cells (%)",
    fill = "Cell subset"
  )


# Save the stacked bar plot.


ggsave(
  "Azimuth_epithelial_ann_level3_stacked_bar.png",
  p_level3_bar,
  width = 9,
  height = 10,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 7. Reference epithelial UMAP
# ============================================================
# Calculate median coordinates for cell-type labels.
ref_epi_centers <- ref_epi_df %>%
  filter(!is.na(ann_level_3)) %>%
  group_by(ann_level_3) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  )

p_ref_epi_alone <- ggplot(
  ref_epi_df,
  aes(x = UMAP_1, y = UMAP_2, color = ann_level_3)
) +
  geom_point(size = 0.18, alpha = 0.65) +
  ggrepel::geom_text_repel(
    data = ref_epi_centers,
    aes(x = UMAP_1, y = UMAP_2, label = ann_level_3),
    inherit.aes = FALSE,
    size = 10,  # Cell-type label size
    fontface = "bold",
    color = "black",
    max.overlaps = Inf,
    seed = 123
  ) +
  coord_fixed(
    xlim = range(ref_epi_df$UMAP_1, na.rm = TRUE),
    ylim = range(ref_epi_df$UMAP_2, na.rm = TRUE)
  ) +
  scale_color_manual(
    values = ref_colors,
    drop = FALSE
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.7),
    
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5  
    ),
    
    axis.title.x = element_text(
      size = 30,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 30,
      face = "bold"
    )
  ) +
  labs(
    title = "HLCA Reference",
    x = "UMAP1",
    y = "UMAP2"
  )

ggsave(
  "HLCA_epithelial_reference_alone.png",
  p_ref_epi_alone,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 8. Query cells colored by passage over the epithelial reference
# ============================================================

p_query_passage_overlay <- ggplot() +
  
  # Reference epithelial cells
  geom_point(
    data = ref_epi_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    color = "grey85",
    size = 0.15,
    alpha = 0.45
  ) +
  
  # Query cells colored by passage
  geom_point(
    data = query_epi_df,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = passage
    ),
    size = 0.5,
    alpha = 0.8
  ) +
  
  # Same coordinate range as reference UMAP
  coord_fixed(
    xlim = range(ref_epi_df$UMAP_1, na.rm = TRUE),
    ylim = range(ref_epi_df$UMAP_2, na.rm = TRUE)
  ) +
  
  theme_classic() +
  
  theme(
    # ========================================================
    # Axis
    # ========================================================
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.7
    ),
    
    # ========================================================
    # Plot title
    # ========================================================
    plot.title = element_text(
      size = 25,
      face = "bold",
      hjust = 0.5
    ),
    
    # ========================================================
    # Axis titles
    # ========================================================
    axis.title.x = element_text(
      size = 30,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 30,
      face = "bold"
    ),
    
    # ========================================================
    # Legend
    # ========================================================
    legend.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 18
    )
  ) +
  
  labs(
    title = "Query",
    x = "UMAP1",
    y = "UMAP2",
    color = "Passage"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  )


# Save the plot.

ggsave(
  "HLCA_epithelial_query_overlay_by_passage.png",
  p_query_passage_overlay,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 9. Query cells by predicted epithelial subtype and passage
# ============================================================
# Prepare query data for plotting.


query_epi_df$predicted_ann_level_3 <- factor(
  query_epi_df$predicted_ann_level_3,
  levels = ref_levels
)

query_epi_df$passage <- factor(
  query_epi_df$passage,
  levels = passage_levels
)

# Remove cells without a level-3 annotation.
query_epi_df_plot <- query_epi_df %>%
  filter(!is.na(predicted_ann_level_3))


# Add dummy rows to retain all reference cell types in the legend.

legend_dummy <- data.frame(
  UMAP_1 = NA_real_,
  UMAP_2 = NA_real_,
  predicted_ann_level_3 = factor(
    ref_levels,
    levels = ref_levels
  )
)


# Plot projected query cells.

p_query_passage_split <- ggplot() +
  
  # HLCA epithelial reference background
  geom_point(
    data = ref_epi_df,
    aes(
      x = UMAP_1,
      y = UMAP_2
    ),
    color = "grey88",
    size = 0.12,
    alpha = 0.4
  ) +
  
  # Query cells colored by predicted epithelial subset
  geom_point(
    data = query_epi_df_plot,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = predicted_ann_level_3
    ),
    size = 0.45,
    alpha = 0.8
  ) +
  
  # Dummy layer to retain all reference cell subsets in legend
  geom_point(
    data = legend_dummy,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      color = predicted_ann_level_3
    ),
    size = 0,
    alpha = 0,
    show.legend = TRUE,
    inherit.aes = FALSE
  ) +
  
  # One panel per passage
  facet_wrap(
    ~passage,
    nrow = 1,
    drop = FALSE,
    labeller = as_labeller(
      c(
        P00 = "Passage 0",
        P05 = "Passage 5",
        P12 = "Passage 12"
      )
    )
  ) +
  
  # Same coordinates as reference UMAP
  coord_fixed(
    xlim = range(ref_epi_df$UMAP_1, na.rm = TRUE),
    ylim = range(ref_epi_df$UMAP_2, na.rm = TRUE)
  ) +
  
  # Same colors as reference
  scale_color_manual(
    values = ref_colors,
    limits = ref_levels,
    breaks = ref_levels,
    drop = FALSE,
    name = "Cell subset"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    strip.background = element_blank(),
    
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.7
    ),
    
    axis.title.x = element_text(
      size = 30,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 30,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = 22,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = 20,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    )
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(
        size = 5,
        alpha = 1
      )
    )
  ) +
  
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )


# Save

ggsave(
  "HLCA_epithelial_query_split_by_passage.png",
  p_query_passage_split,
  width = 15,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
# ============================================================
# 10. Normalize RNA for expression overlays
# ============================================================
DefaultAssay(query_az_epi) <- "RNA"

query_az_epi <- NormalizeData(
  query_az_epi,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ============================================================
# 11. SFTPC, SFTPB, AGER, KRT5 expression by passage
# ============================================================

marker_genes <- c(
  "SFTPC",
  "SFTPB",
  "AGER",
  "KRT5"
)
# ============================================================
# Check marker-gene availability.

genes_present <- intersect(
  marker_genes,
  rownames(query_az_epi)
)


if (length(genes_present) == 0) {
  stop(
    "None of the requested marker genes are present in the query object."
  )
}


# ============================================================
# Fetch normalized expression data.


expression_df <- FetchData(
  query_az_epi,
  vars = genes_present,
  layer = "data"
)

expression_df$cell <- rownames(expression_df)


# ============================================================
# Convert expression data to long format.


expression_long <- expression_df %>%
  pivot_longer(
    cols = all_of(genes_present),
    names_to = "gene",
    values_to = "expression"
  )


# ============================================================
# Add projected UMAP coordinates and passage metadata.


query_expression_df <- expression_long %>%
  left_join(
    query_epi_df %>%
      select(
        cell,
        UMAP_1,
        UMAP_2,
        passage
      ),
    by = "cell"
  )


# ============================================================
# Set gene and passage order.


query_expression_df$gene <- factor(
  query_expression_df$gene,
  levels = marker_genes
)

query_expression_df$passage <- factor(
  query_expression_df$passage,
  levels = passage_levels
)


# ============================================================
# Plot projected query cells.


p_marker_expression_split <- ggplot() +
  
  # ----------------------------------------------------------
# HLCA epithelial reference background
# ----------------------------------------------------------
geom_point(
  data = ref_epi_df,
  aes(
    x = UMAP_1,
    y = UMAP_2
  ),
  color = "grey90",
  size = 0.10,
  alpha = 0.35
) +
  
  # ----------------------------------------------------------
# Query cells colored by gene expression
# ----------------------------------------------------------
geom_point(
  data = query_expression_df,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = expression
  ),
  size = 0.45,
  alpha = 0.85
) +
  
  # ----------------------------------------------------------
# Gene × passage
# ----------------------------------------------------------
facet_grid(
  rows = vars(gene),
  cols = vars(passage),
  drop = FALSE,
  labeller = labeller(
    passage = as_labeller(
      c(
        P00 = "Passage 0",
        P05 = "Passage 5",
        P12 = "Passage 12"
      )
    )
  )
) +
  
  # ----------------------------------------------------------
# Expression color scale
# ----------------------------------------------------------
ggplot2::scale_color_gradientn(
  colours = c(
    "#FEF9E7",
    "#FEC44F",
    "#F16913",
    "#D7301F",
    "#7F0000"
  )
) +
  
  # ----------------------------------------------------------
# Same UMAP coordinates as reference
# ----------------------------------------------------------
coord_fixed(
  xlim = range(
    ref_epi_df$UMAP_1,
    na.rm = TRUE
  ),
  ylim = range(
    ref_epi_df$UMAP_2,
    na.rm = TRUE
  )
) +
  
  theme_classic() +
  
  theme(
    # --------------------------------------------------------
    # Axis
    # --------------------------------------------------------
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    
    axis.title.x = element_text(
      size = 30,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 30,
      face = "bold"
    ),
    
    # --------------------------------------------------------
    # Border around each UMAP
    # --------------------------------------------------------
    panel.border = element_rect(
      color = "grey",
      fill = NA,
      linewidth = 0.7
    ),
    
    # --------------------------------------------------------
    # Remove facet title boxes
    # --------------------------------------------------------
    strip.background = element_blank(),
    
    # Passage / gene labels
    strip.text = element_text(
      size = 22,
      face = "bold"
    ),
    
    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 16
    )
  ) +
  
  labs(
    x = "UMAP1",
    y = "UMAP2"
  )


# Save

ggsave(
  "HLCA_epithelial_marker_by_passage.png",
  p_marker_expression_split,
  width = 15,
  height = 16,
  dpi = 300,
  bg = "white"
)

ggsave(
  "HLCA_epithelial_marker_by_passage.svg",
  p_marker_expression_split,
  width = 15,
  height = 16,
  dpi = 300,
  bg = "white"
)
# ============================================================
# 12. Save final objects
# ============================================================
saveRDS(
  query_az_hlca,
  "query_Azimuth_HLCA_all_cells.rds"
)

saveRDS(
  query_az_epi,
  "query_Azimuth_HLCA_epithelial_level1_normalized.rds"
)


message("Analysis completed.")

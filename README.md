# Lung_AT2_Identification

R workflows for studying human alveolar type 2 (AT2) cell identity, culture, and differentiation, with supporting endothelial cell and fibroblast expression analyses.

## Folder setup

Create below directories locally and add your input data. Datasets should be downloaded seperately.

```text
Lung_AT2_Identification/
├── src/
└── analysis/
    ├── AT2/
    │   ├── SB001_cellranger_count_outs/filtered_feature_bc_matrix/
    │   ├── SB002_cellranger_count_outs/filtered_feature_bc_matrix/
    │   ├── SB003_cellranger_count_outs/filtered_feature_bc_matrix/
    │   └── Azimuth/
    ├── Diff/
    │   └── SB2124_cellranger_multi_outs/per_sample_outs/
    ├── GSE135893/
    ├── GSE136831/
    └── bulkRNAseq/
```

- **Cell Ranger:** keep the original matrix files. Under `Diff/SB2124_cellranger_multi_outs/per_sample_outs/`, place each sample's matrix at `<sample>/count/sample_filtered_feature_bc_matrix/` for `SB021`, `SB022`, `SB023`, and `SB024`.
- **GEO references:** put each dataset's count matrix, gene list, barcodes, and metadata in its respective `GSE135893/` or `GSE136831/` folder, retaining the filenames specified in scripts `03` and `04`.
- **Bulk RNA-seq:** put `Set1_pooled_edgeR_analysis.RData` (containing the raw-count `dge_set1` object) and `metadata_set1(1)(1).csv` (with `sample` and `condition` columns) in `bulkRNAseq/`.

### Acknowledgement

Part of scripts was generated using generative coding agent.
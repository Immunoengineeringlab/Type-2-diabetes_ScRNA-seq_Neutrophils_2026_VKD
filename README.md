# scRNA-seq of Peripheral Blood Leukocytes: Type 2 Diabetes vs Healthy

**ArrayExpress / BioStudies accession:** E-MTAB-XXXXX  
**Deposited by:** Vinod Kumar Dorai, Jhunjhunwala Lab, Department of Bioengineering, Indian Institute of Science (IISc), Bengaluru, India  
**Contact:** siddharth@iisc.ac.in  
**Associated publication:** Dorai VK et al. *in review* [year]. DOI: [XXXX]  

---

## Important Note on Dataset Scope

This sequencing run was part of a larger multiplexed study that also included participants with **chronic kidney disease**. The data for those participants belong to a separate, ongoing study and will be deposited independently upon publication of that work.

The present deposit contains the **complete raw sequencing and processed data** as generated from the BD Rhapsody pipeline for the entire multiplexed run.

**To reproduce the analyses described in the associated publication, restrict analysis to the following four sample tags only:**

| Sample Tag | Condition | Sex | Age (yr) | Cartridge |
|---|---|---|---|---|
| SampleTag01_hs | Type 2 diabetes (T2D) | Male | 57 | Cart3 |
| SampleTag04_hs | Healthy | Male | 59 | Cart3 |
| SampleTag05_hs | Type 2 diabetes (T2D) | Female | 56 | Cart3 |
| SampleTag12_hs | Healthy | Female | 52 | Cart4 |

All other sample tags correspond to the separate ongoing study. The analysis scripts implement this subsetting explicitly via the `VALID_TAGS` vector in `01_load_data_1.R`.

---

## Study Overview

Single-cell RNA sequencing of peripheral blood leukocytes from apparently healthy and type 2 diabetic (T2D) donors, generated to characterise whether the circulating immune compartment is transcriptionally altered in T2D. The study identified **neutrophils as the most transcriptionally remodelled immune population**, with 618 differentially expressed genes (adjusted P < 0.05, |log2FC| > 0.25), and confirmed that the core neutrophil signature was reproducible across individual donors.

| Parameter | Detail |
|---|---|
| Platform | BD Rhapsody (whole-transcriptome WTA + AbSeq antibody-derived tags) |
| Cartridges | 2 (Cart3, Cart4), sample-tag multiplexed |
| Sequencing | Illumina Novaseq 6000 platform, Eurofins Genomics, Bengaluru, India |
| Raw processing | BD Rhapsody RSEC pipeline on Seven Bridges Genomics |
| Downstream analysis | R (>= 4.4.0), Seurat v5.1.0 |

---

## Deposited Files

### Raw Data (complete cartridge output)

> Raw data contains all multiplexed samples.

### Processed Data (4 T2D/healthy donors only)

```
CART3HUMAN_RSEC_MolsPerCell_MEX/
  cart3_matrix.mtx.gz        Gene-expression count matrix (RSEC-adjusted, Cart3)
  cart3_features.tsv.gz      Gene names
  cart3_barcodes.tsv.gz      Cell barcodes
CART3HUMAN_Sample_Tag_Calls.csv   Per-cell sample-tag assignments, Cart3

CART4HUMAN_RSEC_MolsPerCell_MEX/
  cart4_matrix.mtx.gz        Gene-expression count matrix (RSEC-adjusted, Cart4)
  cart4_features.tsv.gz      Gene names
  cart4_barcodes.tsv.gz      Cell barcodes
CART4HUMAN_Sample_Tag_Calls.csv   Per-cell sample-tag assignments, Cart4
```

### Analysis Scripts

| Script | Description |
|---|---|
| `01_load_data_1.R` | Load raw MEX + tag CSVs; subset to 4 donors |
| `02_qc_filter_1.R` | QC metrics and cell filtering |
| `03_normalise_pca_1.R` | SCTransform normalisation and PCA |
| `04_harmony_cluster_1.R` | Harmony batch correction, clustering, UMAP |
| `05_annotate_singler_1.R` | Cell-type annotation (SingleR, Monaco reference) |
| `06_cell_proportions_1.R` | Cell-type proportion visualisation (scMEGA) |
| `07_deg_analysis_1.R` | Differential expression per cell type (MAST) |
| `07b_neutrophil_DEG_donor_consistency.R` | Per-donor reproducibility check for neutrophil DEGs |
| `08_deg_visualisation_1.R` | DEG figures including reproducible-gene volcano |
| `09_deg_downstream.R` | Pathway enrichment (GO, KEGG, GSEA) |

---

## Reproducing the Published Analysis

Scripts are numbered in execution order (`01` → `09`). Run them sequentially in the same working directory. Each script loads the checkpoint saved by the previous one.

### Prerequisites

```r
# Bioconductor packages
BiocManager::install(c("SingleR", "celldex", "clusterProfiler",
                       "org.Hs.eg.db", "enrichplot", "ReactomePA"))

# CRAN packages
install.packages(c("Seurat", "harmony", "ggplot2", "dplyr", "tidyr",
                   "patchwork", "ggrepel", "RColorBrewer", "pheatmap",
                   "scales", "stringr", "viridis", "ggtext", "tibble",
                   "msigdbr", "fgsea", "Matrix", "scMEGA"))
```

`set.seed(42)` is used throughout for reproducibility.

---

### Section 01 — Load Data (`01_load_data_1.R`)

Loads Cart3 and Cart4 MEX directories and `Sample_Tag_Calls` CSVs. Creates a merged Seurat object with Gene Expression (RNA assay) and Antibody Capture (Protein assay). Barcodes are prefixed `C3_` or `C4_` to prevent collision between cartridges. Subsets to the four T2D/healthy sample tags and assigns Condition metadata.

**Update the path variables at the top of the script to match your local file locations.**

- Saves: `Merged.rds`
- Output cells: **7,906** (SampleTag01=2108, SampleTag04=1765, SampleTag05=1350, SampleTag12=2683)

---

### Section 02 — QC Filtering (`02_qc_filter_1.R`)

Computes per-cell QC metrics (mitochondrial, ribosomal, haemoglobin, heat-shock, platelet gene content) and applies hard filters:

| Filter | Threshold |
|---|---|
| nFeature_RNA | 300 – 4,000 genes per cell |
| nCount_RNA | 500 – 20,000 transcripts per cell |
| percent.mt | < 20% |
| percent.hb | < 1% |

- Saves: `Merged_after_qc.rds`

---

### Section 03 — Normalisation + PCA (`03_normalise_pca_1.R`)

`JoinLayers` is called before SCTransform (required for Seurat v5 merged objects). SCTransform v2 (`vst.flavor = "v2"`) is applied, regressing `percent.mt` only. Ribosomal content is not regressed as it is biologically meaningful in immune cells. Cartridge batch is corrected downstream by Harmony, not here.

PC selection is empirical: minimum of:
- `co1`: first PC where cumulative variance > 90% and individual variance < 5%
- `co2`: last PC where change in variance between consecutive PCs > 0.1%

---

### Section 04 — Harmony Integration + Clustering (`04_harmony_cluster_1.R`)

Harmony corrects for both `Cartridge` (technical batch) and `Condition` (preserving biological signal). Clustering uses Louvain algorithm. Cartridge mixing is verified per cluster as a quality check.

| Parameter | Value |
|---|---|
| group.by.vars | Cartridge + Condition |
| k.param | 20 |
| resolution | 0.5 |

- Saves: `merged_harmony.rds`

---

### Section 05 — Cell-Type Annotation (`05_annotate_singler_1.R`)

The RNA assay is log-normalised (`NormalizeData`) before SingleR, which requires log-normalised input. SingleR is run against `celldex::MonacoImmuneData()` for both main and fine labels. Pruned labels are used, setting low-confidence calls to `NA`. Both `Cell_annot` (main) and `Cell_annot_fine` columns are added to metadata.

- Saves: `merged_annotated.rds`

---

### Section 06 — Cell Proportions (`06_cell_proportions_1.R`)

Uses `CellPropPlot` from the scMEGA package to visualise and export cell-type proportions per condition.

- Saves: `cell_type_proportions.csv`

---

### Section 07 — Differential Expression (`07_deg_analysis_1.R`)

`PrepSCTFindMarkers` recalculates corrected SCT counts prior to `FindMarkers`. DEG testing is performed per cell type.

| Parameter | Value |
|---|---|
| test.use | MAST |
| assay | SCT |
| min.pct | 0.10 |
| ident.1 / ident.2 | Diabetic / Healthy |
| Adjusted P threshold | < 0.05 (Benjamini-Hochberg) |
| \|log2FC\| threshold | > 0.25 |
| Minimum cells | ≥ 20 total; ≥ 10 per condition |

- Saves: `merged_with_DEGs.rds`, `DEG_summary_all_celltypes.csv`, `DEG_full_results_all_celltypes.csv`

---

### Section 07b — Per-Donor Reproducibility (`07b_neutrophil_DEG_donor_consistency.R`)

Tests whether the cohort-level neutrophil DEG signature is reproducible across donors rather than driven by a single individual. For each of the four cross-donor pairings (each T2D donor vs each healthy donor), DEG testing is re-run with the same MAST parameters as Section 07.

For each cohort-level DEG, consistency is scored as:
- **Strict:** same direction AND adjusted P < 0.05 in each pairing
- **Lenient:** same direction of log2FC regardless of per-pairing significance

The lenient metric is recommended for interpretation given the reduced per-pairing statistical power relative to the full cohort comparison. This analysis generates the per-donor heatmap (Fig. 1E in the manuscript). All figures use `pheatmap(filename=)` or `ggsave()` — no `dev.off()` calls.

- Saves: `neutrophil_DEG_donor_consistency.csv`, `neutrophil_donor_mean_expression.csv`, `neutrophil_DEG_per_pairing.csv`, `neutrophil_DEG_per_pairing_counts.csv`, `Fig_neutrophil_DEG_consistency_heatmap.pdf`, `Fig_neutrophil_donor_pairing_concordance.pdf`

---

### Section 08 — DEG Visualisation (`08_deg_visualisation_1.R`)

Generates publication figures from the DEG results CSVs. The neutrophil volcano plot labels only genes confirmed as reproducible across donors in Section 07b (`REPRODUCIBLE_GENES` vector). Several interferon-stimulated genes showed large cohort-level fold changes but were concentrated in a single donor and are therefore not highlighted. All DEG points are plotted; reproducible genes are visually emphasised.

---

### Section 09 — Pathway Enrichment (`09_deg_downstream.R`)

GO (biological process) and KEGG over-representation analysis via `clusterProfiler`. GSEA via `fgsea` using GO:BP gene sets from `msigdbr` (Homo sapiens, C5/GO:BP). Enrichment is run on the full neutrophil DEG list ranked by log2FC.

> **Note:** `org.Hs.eg.db` and `AnnotationDbi` mask `dplyr` verbs (`select`, `filter`, `rename`, `mutate`). These are explicitly reclaimed using `dplyr::` prefixes at the top of the script.

---

## Key Parameter Summary

| Step | Parameter | Value |
|---|---|---|
| Reproducibility | set.seed | 42 |
| QC | nFeature_RNA | 300 – 4,000 |
| QC | nCount_RNA | 500 – 20,000 |
| QC | percent.mt | < 20% |
| QC | percent.hb | < 1% |
| Normalisation | Method | SCTransform v2 |
| Normalisation | Regression | percent.mt only |
| Integration | Method | Harmony |
| Integration | group.by.vars | Cartridge + Condition |
| Clustering | Algorithm | Louvain |
| Clustering | Resolution | 0.5 |
| Clustering | k | 20 |
| Annotation | Reference | Monaco Immune Data (celldex) |
| Annotation | Labels | pruned.labels (main + fine) |
| DEG test | Method | MAST |
| DEG test | Assay | SCT |
| DEG test | min.pct | 0.10 |
| DEG thresholds | Adjusted P | < 0.05 (BH) |
| DEG thresholds | \|log2FC\| | > 0.25 |
| Enrichment | Tools | clusterProfiler, fgsea |
| Enrichment | Gene sets | GO:BP (msigdbr C5/GO:BP), KEGG |

---

## Licence

Data: [Creative Commons Attribution 4.0 (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)  
Scripts: [MIT licence](https://opensource.org/licenses/MIT)

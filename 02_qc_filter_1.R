# ============================================================
# SECTION 02: QC METRICS & CELL FILTERING
# ============================================================
# INPUT  : merged  (from Section 01, or reload checkpoint below)
# OUTPUT : merged  (QC-filtered Seurat object)
# SAVES  : Merged_after_qc.rds
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# ── Reload checkpoint if running standalone ────────────────────────────────
merged <- readRDS("Merged.rds")   # skip to next section instead

# ============================================================
# 2.1  COMPUTE QC METRICS
# ============================================================

DefaultAssay(merged) <- "RNA"

merged[["percent.mt"]]   <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")    # fixed: was ^RPS|^RPL
merged[["percent.hb"]]   <- PercentageFeatureSet(merged, pattern = "^HB[^BP]")   # fixed: was ^HB[^(P)]
merged[["percent.hsp"]]  <- PercentageFeatureSet(merged, pattern = "^HSP|^DNAJ")
merged[["percent.plat"]] <- PercentageFeatureSet(merged, pattern = "PECAM1|PF4")

cat("=== QC Summary (pre-filter) ===\n")
print(summary(merged@meta.data[, c("nFeature_RNA", "nCount_RNA",
                                   "percent.mt", "percent.ribo",
                                   "percent.hb")]))

# ============================================================
# 2.2  QC VISUALISATIONS (pre-filter)
# ============================================================

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt",
                 "percent.ribo", "percent.hb",
                 "percent.hsp",  "percent.plat")

# Per sample tag
VlnPlot(merged,
        features = qc_features,
        group.by = "Sample_Tag",
        pt.size  = 0,
        ncol     = 4) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.title  = element_text(face = "bold", size = 9))

# Per cartridge — detect batch-level differences
VlnPlot(merged,
        features = c("nFeature_RNA", "nCount_RNA",
                     "percent.mt",   "percent.ribo"),
        group.by = "Cartridge",
        pt.size  = 0,
        cols     = c("Cart3" = "#1D6FA4", "Cart4" = "#E36B2D"),
        ncol     = 4) +
  theme(plot.title = element_text(face = "bold", size = 9)) +
  plot_annotation(title = "QC Metrics per Cartridge (pre-filter)")

# Scatter plots
p_scatter1 <- FeatureScatter(merged,
                             feature1 = "nCount_RNA",
                             feature2 = "percent.mt",
                             group.by = "Cartridge") +
  ggtitle("Count vs MT%")

p_scatter2 <- FeatureScatter(merged,
                             feature1 = "nCount_RNA",
                             feature2 = "nFeature_RNA",
                             group.by = "Cartridge") +
  ggtitle("Count vs Features")

print(p_scatter1 | p_scatter2)

# ============================================================
# 2.3  APPLY QC FILTERS
# ============================================================

cat("\n=== Applying QC Filters ===\n")
cat("nFeature_RNA : 300 – 4000\n")
cat("nCount_RNA   : 500 – 20000\n")
cat("percent.mt   : < 20\n")
cat("percent.hb   : < 1\n\n")

cells_before <- ncol(merged)

merged <- subset(
  merged,
  subset = nFeature_RNA > 300  &
           nFeature_RNA < 4000 &
           nCount_RNA   > 500  &
           nCount_RNA   < 20000 &
           percent.mt   < 20   &
           percent.hb   < 1
)

cells_after <- ncol(merged)

cat("Cells before filter :", cells_before, "\n")
cat("Cells after  filter :", cells_after,  "\n")
cat("Cells removed       :", cells_before - cells_after, "\n")

cat("\n=== Post-filter: Cells per Cartridge ===\n")
print(table(merged$Cartridge))
cat("\n=== Post-filter: Cells per Sample_Tag ===\n")
print(table(merged$Sample_Tag))
cat("\n=== Post-filter: Cells per Condition ===\n")
print(table(merged$Condition))

# ============================================================
# 2.4  POST-FILTER VISUALISATIONS
# ============================================================

VlnPlot(merged,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "Cartridge",
        pt.size  = 0,
        cols     = c("Cart3" = "#1D6FA4", "Cart4" = "#E36B2D"),
        ncol     = 3) +
  plot_annotation(title = "QC Metrics per Cartridge (post-filter)")

print(
  FeatureScatter(merged,
                 feature1 = "nCount_RNA",
                 feature2 = "nFeature_RNA",
                 group.by = "Condition") +
    ggtitle("Count vs Features — post-filter")
)

print(summary(merged@meta.data[, c("nFeature_RNA", "nCount_RNA",
                                   "percent.mt")]))

# ============================================================
# 2.5  SAVE CHECKPOINT
# ============================================================

saveRDS(merged, "Merged_after_qc.rds")
cat("\n✅ Saved: Merged_after_qc.rds\n")
cat("✅ Section 02 complete\n")

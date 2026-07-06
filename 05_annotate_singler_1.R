# ============================================================
# SECTION 05: CELL TYPE ANNOTATION — SingleR
# ============================================================
# INPUT  : merged_harmony.rds
# OUTPUT : merged  (with Cell_annot and Cell_annot_fine columns)
#          propcolors  (named colour vector for downstream plots)
# SAVES  : merged_annotated.rds
# ============================================================

library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(RColorBrewer)

# ── Load checkpoint ────────────────────────────────────────────────────────
merged <- readRDS("merged_harmony.rds")
cat("Loaded: merged_harmony.rds |", ncol(merged), "cells\n")

# ============================================================
# 5.1  PREPARE LOG-NORMALISED RNA FOR SingleR
# ============================================================

DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)
merged <- NormalizeData(merged, assay = "RNA", verbose = FALSE)

norm_counts <- GetAssayData(merged, assay = "RNA", layer = "data")

cat("Expression matrix for SingleR:", dim(norm_counts), "\n")

# ============================================================
# 5.2  LOAD REFERENCE & RUN SingleR
# Monaco Immune reference — appropriate for human PBMCs
# ============================================================

ref <- celldex::MonacoImmuneData()

cat("Running SingleR (main labels)...\n")
results_main <- SingleR(
  test   = norm_counts,
  ref    = ref,
  labels = ref$label.main
)

cat("Running SingleR (fine labels)...\n")
results_fine <- SingleR(
  test   = norm_counts,
  ref    = ref,
  labels = ref$label.fine
)

# ============================================================
# 5.3  ADD ANNOTATIONS TO METADATA
# Using pruned.labels — removes low-confidence ambiguous calls
# ============================================================

merged <- AddMetaData(merged,
                      results_main$pruned.labels,
                      col.name = "Cell_annot")
merged <- AddMetaData(merged,
                      results_fine$pruned.labels,
                      col.name = "Cell_annot_fine")

cat("\n=== Annotation Summary (main) ===\n")
print(table(merged$Cell_annot, useNA = "ifany"))

cat("\n=== Ambiguous cells (pruned to NA) ===\n")
cat("Main:", sum(is.na(results_main$pruned.labels)),
    "/", ncol(merged), "\n")
cat("Fine:", sum(is.na(results_fine$pruned.labels)),
    "/", ncol(merged), "\n")

# ============================================================
# 5.4  ANNOTATION QUALITY DIAGNOSTICS
# ============================================================

print(
  plotScoreHeatmap(results_main) +
    ggtitle("SingleR Score Heatmap — Main Labels")
)

print(
  plotDeltaDistribution(results_main, ncol = 4) +
    ggtitle("Delta Distribution — Low scores = ambiguous")
)

# ============================================================
# 5.5  ANNOTATED UMAPs
# ============================================================

nb_cols    <- length(unique(na.omit(merged$Cell_annot)))
propcolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb_cols)
names(propcolors) <- sort(unique(na.omit(merged$Cell_annot)))

p_annot_main <- DimPlot(merged,
                        reduction = "umap_integrated",
                        group.by  = "Cell_annot",
                        label     = TRUE, repel = TRUE,
                        pt.size   = 0.6,
                        cols      = propcolors) +
  ggtitle("Cell Type Annotation — Main (SingleR)")

p_annot_fine <- DimPlot(merged,
                        reduction = "umap_integrated",
                        group.by  = "Cell_annot_fine",
                        label     = TRUE, repel = TRUE,
                        pt.size   = 0.6) +
  ggtitle("Cell Type Annotation — Fine (SingleR)")

print(p_annot_main)
print(p_annot_fine)


Umap_annot <-   DimPlot(merged,
          reduction = "umap_integrated",
          group.by  = "Cell_annot",
          cols      = propcolors,
          pt.size   = 0.3,
          ncol      = 2) +
    ggtitle("Cell Types — Healthy vs Diabetic")

ggsave("Fig1_UMAP.pdf", Umap_annot,
       width = 10, height = 5.5, dpi = 300)
# ============================================================
# 5.6  SAVE CHECKPOINT
# ============================================================

saveRDS(merged, "merged_annotated.rds")
cat("\n✅ Saved: merged_annotated.rds\n")
cat("✅ Section 05 complete\n")

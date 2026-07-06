# ============================================================
# SECTION 07: DEG ANALYSIS — All Cell Types
# Diabetic vs Healthy | SCT | Wilcoxon
# ============================================================
# INPUT  : merged_annotated.rds
# OUTPUT : summary_df   (DEG counts per cell type)
#          all_degs     (full DEG results table)
# SAVES  : merged_with_DEGs.rds
#          DEG_summary_all_celltypes.csv
#          DEG_full_results_all_celltypes.csv
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
setwd("C:/IISc/ScRNA_Analysis_Human/20260622_output")

# ── Load checkpoint ────────────────────────────────────────────────────────
merged <- readRDS("merged_annotated.rds")
cat("Loaded: merged_annotated.rds |", ncol(merged), "cells\n")

# ============================================================
# 7.1  PREPARE SCT FOR FindMarkers
# PrepSCTFindMarkers recalculates corrected counts across conditions
# Must be run before FindMarkers with assay = "SCT"
# ============================================================

DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)

DefaultAssay(merged) <- "SCT"

cat("Running PrepSCTFindMarkers...\n")
merged <- PrepSCTFindMarkers(merged, verbose = TRUE)
cat("✅ PrepSCTFindMarkers complete\n")

# ============================================================
# 7.2  IDENTIFY ELIGIBLE CELL TYPES
# Minimum 20 cells total; enforced per-condition (10 each) in loop
# ============================================================

Idents(merged) <- "Cell_annot"

cell_types <- names(table(merged$Cell_annot))[
  table(merged$Cell_annot) >= 20
]

cat("\n=== Cell types eligible for DEG analysis ===\n")
print(cell_types)

# Per-condition cell count breakdown
count_df <- merged@meta.data %>%
  filter(!is.na(Cell_annot)) %>%
  group_by(Cell_annot, Condition) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  pivot_wider(names_from  = Condition,
              values_from = n_cells,
              values_fill = 0) %>%
  mutate(Total = Healthy + Diabetic) %>%
  arrange(desc(Total))

cat("\n=== Cell counts per type per condition ===\n")
print(count_df)

# ============================================================
# 7.3  DEG TESTING: Diabetic vs Healthy per cell type
# ============================================================

deg_results <- list()   # full marker tables
deg_summary <- list()   # UP / DOWN / total counts

cat("\n=== Running DEG per cell type ===\n")

for (ct in cell_types) {

  cat(sprintf("  Processing: %-30s", ct))

  ct_cells   <- subset(merged, idents = ct)
  n_healthy  <- sum(ct_cells$Condition == "Healthy")
  n_diabetic <- sum(ct_cells$Condition == "Diabetic")

  if (n_healthy < 10 || n_diabetic < 10) {
    cat(sprintf("SKIPPED (Healthy=%d, Diabetic=%d)\n",
                n_healthy, n_diabetic))
    next
  }

  Idents(ct_cells) <- "Condition"

  tryCatch({
    markers <- FindMarkers(
      ct_cells,
      ident.1         = "Diabetic",
      ident.2         = "Healthy",
      assay           = "SCT",
      test.use        = "MAST",
      min.pct         = 0.10,
      verbose         = FALSE
    )

    markers$gene      <- rownames(markers)
    markers$cell_type <- ct
    rownames(markers) <- NULL

    markers$direction <- dplyr::case_when(
      markers$avg_log2FC >  0.25 & markers$p_val_adj < 0.05 ~ "Upregulated",
      markers$avg_log2FC < -0.25 & markers$p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )

    deg_results[[ct]] <- markers

    n_up    <- sum(markers$direction == "Upregulated")
    n_down  <- sum(markers$direction == "Downregulated")
    n_total <- n_up + n_down

    deg_summary[[ct]] <- data.frame(
      Cell_Type  = ct,
      N_Healthy  = n_healthy,
      N_Diabetic = n_diabetic,
      UP         = n_up,
      DOWN       = n_down,
      Total_DEGs = n_total,
      stringsAsFactors = FALSE
    )

    cat(sprintf("UP=%d | DOWN=%d | TOTAL=%d\n", n_up, n_down, n_total))

  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
  })
}

# ============================================================
# 7.4  COMBINE RESULTS
# ============================================================

summary_df <- dplyr::bind_rows(deg_summary) %>%
  arrange(desc(Total_DEGs))

all_degs <- dplyr::bind_rows(deg_results)

cat("\n==============================\n")
cat("=== DEG SUMMARY TABLE ===\n")
cat("==============================\n")
print(summary_df)

# ============================================================
# 7.5  SAVE OUTPUTS
# ============================================================

saveRDS(merged, "merged_with_DEGs.rds")
write.csv(summary_df, "DEG_summary_all_celltypes.csv",     row.names = FALSE)
write.csv(all_degs,   "DEG_full_results_all_celltypes.csv", row.names = FALSE)

cat("\n✅ Saved: merged_with_DEGs.rds\n")
cat("✅ Saved: DEG_summary_all_celltypes.csv\n")
cat("✅ Saved: DEG_full_results_all_celltypes.csv\n")
cat("✅ Section 07 complete\n")

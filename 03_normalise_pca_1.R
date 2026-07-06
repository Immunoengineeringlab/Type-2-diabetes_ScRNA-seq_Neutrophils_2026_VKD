# ============================================================
# SECTION 03: NORMALISATION (SCTransform) + PCA + PRE-HARMONY UMAP
# ============================================================
# INPUT  : Merged_after_qc.rds
# OUTPUT : merged  (SCT normalised, PCA computed, pre-Harmony UMAP)
#          pcs     (integer — selected PC count, used in Section 04)
# SAVES  : (none — heavy compute; checkpoint saved in Section 04)
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# ── Load checkpoint ────────────────────────────────────────────────────────
merged <- readRDS("Merged_after_qc.rds")
cat("Loaded: Merged_after_qc.rds |", ncol(merged), "cells\n")

# ============================================================
# 3.1  SCTRANSFORM NORMALISATION
# Regression: percent.mt ONLY
# percent.ribo NOT regressed — biologically meaningful in immune cells
# Cartridge batch handled downstream by Harmony — do NOT regress here
# ============================================================

DefaultAssay(merged) <- "RNA"

# orig.ident must reflect Cartridge for SCT model fitting per batch
merged$orig.ident <- factor(as.character(merged$Cartridge))
cat("orig.ident levels:", levels(merged$orig.ident), "\n")

# Join layers required before SCTransform in Seurat v5
merged <- JoinLayers(merged)
cat("RNA layers after join:\n")
print(Layers(merged[["RNA"]]))

merged <- SCTransform(
  merged,
  vst.flavor      = "v2",
  vars.to.regress = "percent.mt",
  verbose         = TRUE,
  seed.use        = 42
)

cat("\n✅ SCTransform complete\n")
cat("SCT models stored:", length(merged[["SCT"]]@SCTModel.list), "\n")

# ============================================================
# 3.2  PCA
# ============================================================

DefaultAssay(merged) <- "SCT"

merged <- RunPCA(
  merged,
  verbose        = FALSE,
  reduction.name = "pca_rna_unintegrated"
)

# ── Empirical PC selection ─────────────────────────────────────────────────
pct  <- merged[["pca_rna_unintegrated"]]@stdev /
        sum(merged[["pca_rna_unintegrated"]]@stdev) * 100
cumu <- cumsum(pct)

co1  <- which(cumu > 90 & pct < 5)[1]
co2  <- sort(which((pct[1:(length(pct) - 1)] -
                     pct[2:length(pct)]) > 0.1),
             decreasing = TRUE)[1] + 1
pcs  <- min(co1, co2)

message("Selected number of PCs: ", pcs)

print(
  ElbowPlot(merged,
            reduction = "pca_rna_unintegrated",
            ndims     = 50) +
    geom_vline(xintercept = pcs,
               linetype = "dashed", color = "red") +
    ggtitle(paste0("Elbow Plot — Selected PCs: ", pcs))
)

# ============================================================
# 3.3  PRE-HARMONY UMAP (batch structure visualisation)
# Used only to confirm batch effect exists before Harmony
# ============================================================

merged <- FindNeighbors(merged,
                        reduction = "pca_rna_unintegrated",
                        dims      = 1:pcs,
                        verbose   = FALSE)
merged <- FindClusters(merged,
                       resolution = 0.5,
                       verbose    = FALSE)
merged <- RunUMAP(merged,
                  reduction      = "pca_rna_unintegrated",
                  dims           = 1:pcs,
                  reduction.name = "umap_unintegrated",
                  verbose        = FALSE)

p_pre_cart <- DimPlot(merged,
                      reduction = "umap_unintegrated",
                      group.by  = "Cartridge",
                      pt.size   = 0.3) +
  ggtitle("PRE-Harmony — Cartridge")

p_pre_cond <- DimPlot(merged,
                      reduction = "umap_unintegrated",
                      group.by  = "Condition",
                      pt.size   = 0.3) +
  ggtitle("PRE-Harmony — Condition")

p_pre_tag  <- DimPlot(merged,
                      reduction = "umap_unintegrated",
                      group.by  = "Sample_Tag",
                      pt.size   = 0.3) +
  ggtitle("PRE-Harmony — Sample Tag")

print(p_pre_cart | p_pre_cond | p_pre_tag)

cat("\n✅ Section 03 complete — 'merged' and 'pcs' ready for Section 04\n")

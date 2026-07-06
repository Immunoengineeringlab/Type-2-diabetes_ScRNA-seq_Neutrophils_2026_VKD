# ============================================================
# SECTION 04: HARMONY INTEGRATION + CLUSTERING + POST-HARMONY UMAP
# ============================================================
# INPUT  : merged  (from Section 03, with pca_rna_unintegrated)
#          pcs     (integer from Section 03)
# OUTPUT : merged  (Harmony-corrected, clustered, UMAP computed)
# SAVES  : merged_harmony.rds
# ============================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(patchwork)

# ── If running standalone, reload and recompute pcs ───────────────────────
# merged <- readRDS("Merged_after_qc.rds")
# Then re-run SCTransform + PCA from Section 03, or save merged post-PCA
# and reload here. Easiest: save merged after PCA with saveRDS().

# ============================================================
# 4.1  HARMONY BATCH CORRECTION
# group.by.vars: Cartridge (technical) + Condition (biological)
# Cartridge removes inter-cartridge batch
# Condition ensures Healthy/Diabetic biology is NOT collapsed
# ============================================================

# Ensure factor levels are clean before Harmony
merged$Cartridge  <- factor(as.character(merged$Cartridge))
merged$Condition  <- factor(as.character(merged$Condition))
merged$orig.ident <- factor(as.character(merged$Cartridge))

cat("Cartridge levels:", levels(merged$Cartridge), "\n")
cat("Condition levels:", levels(merged$Condition),  "\n")

merged <- RunHarmony(
  object           = merged,
  group.by.vars    = c("Cartridge", "Condition"),
  reduction        = "pca_rna_unintegrated",
  reduction.save   = "pca_rna_integrated",
  dims.use         = 1:pcs,
  plot_convergence = TRUE,
  verbose          = TRUE
)

cat("✅ Harmony integration complete\n")

# ============================================================
# 4.2  POST-HARMONY CLUSTERING + UMAP
# ============================================================

merged <- FindNeighbors(merged,
                        reduction = "pca_rna_integrated",
                        dims      = 1:pcs,
                        k.param   = 20,
                        verbose   = FALSE)
merged <- FindClusters(merged,
                       resolution = 0.5,
                       verbose    = FALSE)
merged <- RunUMAP(merged,
                  reduction      = "pca_rna_integrated",
                  dims           = 1:pcs,
                  reduction.name = "umap_integrated",
                  verbose        = FALSE)

cat("Clusters found:", nlevels(merged$seurat_clusters), "\n")

# ============================================================
# 4.3  VISUALISATIONS — Before vs After Harmony
# ============================================================

p_post_cart <- DimPlot(merged,
                       reduction = "umap_integrated",
                       group.by  = "Cartridge",
                       pt.size   = 0.3) +
  ggtitle("POST-Harmony — Cartridge")

p_post_cond <- DimPlot(merged,
                       reduction = "umap_integrated",
                       group.by  = "Condition",
                       pt.size   = 0.3) +
  ggtitle("POST-Harmony — Condition")

p_post_clust <- DimPlot(merged,
                        reduction  = "umap_integrated",
                        group.by   = "seurat_clusters",
                        label      = TRUE,
                        label.size = 4,
                        repel      = TRUE,
                        pt.size    = 0.3) +
  ggtitle("POST-Harmony — Clusters")

# Side-by-side batch correction validation
print(
  (p_pre_cart | p_post_cart) /
  (p_pre_cond | p_post_cond)
)
print(p_post_clust)

print(
  DimPlot(merged,
          reduction = "umap_integrated",
          group.by  = "seurat_clusters",
          split.by  = "Condition",
          label     = TRUE,
          pt.size   = 0.3,
          ncol      = 2) +
    ggtitle("Clusters by Condition — Post-Harmony")
)

# ============================================================
# 4.4  CARTRIDGE MIXING QC PER CLUSTER
# Good integration: each cluster ~50/50 from both cartridges
# ============================================================

mix_df <- merged@meta.data %>%
  group_by(seurat_clusters, Cartridge) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

print(
  ggplot(mix_df,
         aes(x = seurat_clusters, y = pct, fill = Cartridge)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 50,
               linetype = "dashed", color = "white", linewidth = 0.8) +
    scale_fill_manual(values = c("Cart3" = "#1D6FA4",
                                 "Cart4" = "#E36B2D")) +
    labs(title    = "Cartridge Mixing per Cluster (Post-Harmony)",
         subtitle = "Good integration ≈ 50/50 per cluster",
         x = "Cluster", y = "Proportion (%)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
)

poor_mix <- mix_df %>% filter(pct > 80)
if (nrow(poor_mix) > 0) {
  cat("⚠️  Poorly mixed clusters (>80% one cartridge):\n")
  print(poor_mix)
} else {
  cat("✅ All clusters well mixed across cartridges\n")
}

# ============================================================
# 4.5  SAVE CHECKPOINT
# ============================================================

saveRDS(merged, "merged_harmony.rds")
cat("\n✅ Saved: merged_harmony.rds\n")
cat("✅ Section 04 complete\n")

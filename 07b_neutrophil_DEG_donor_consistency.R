# ============================================================
# SECTION 07b: PER-DONOR CONSISTENCY OF NEUTROPHIL DEGs
# Purpose: test whether the cohort-level neutrophil signature
#          (Diabetic vs Healthy) is reproducible across donors,
#          or driven by a single outlier donor.
#
# Design note: with 2 healthy + 2 diabetic donors, this is a
#   REPRODUCIBILITY / CONSISTENCY check, NOT a group-level
#   significance test. We do not compute a between-donor p-value.
#   We (i) re-derive DE for each of the 4 cross-donor pairings,
#   (ii) measure how often each cohort-level DEG is recovered in
#   the same direction, and (iii) check per-donor mean expression
#
# INPUT  : merged_with_DEGs.rds        (from 07)
#          DEG_full_results_all_celltypes.csv   (from 07; for the
#                                        cohort-level neutrophil DEG list)
# OUTPUT : neutrophil_DEG_donor_consistency.csv
#          neutrophil_donor_mean_expression.csv
#          Fig_neutrophil_DEG_consistency_heatmap.pdf
#          Fig_neutrophil_donor_pairing_concordance.pdf
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# ── 0. CONFIG ───────────────────────────────────────────────
NEUTROPHIL_LABEL <- "Neutrophils"   # <-- set to the exact Cell_annot value
#     for neutrophils in your object
LFC_CUT  <- 0.25                    # same thresholds as script 07
PADJ_CUT <- 0.05
MIN_PCT  <- 0.10
TEST_USE <- "MAST"

# Donor (sample-tag) -> condition map, from 01_load_data
DIABETIC_TAGS <- c("SampleTag01_hs", "SampleTag05_hs")
HEALTHY_TAGS  <- c("SampleTag04_hs", "SampleTag12_hs")

# ── 1. LOAD ─────────────────────────────────────────────────
merged <- readRDS("merged_with_DEGs.rds")
cat("Loaded merged_with_DEGs.rds |", ncol(merged), "cells\n")

# Sanity: confirm the neutrophil label exists
if (!NEUTROPHIL_LABEL %in% unique(merged$Cell_annot)) {
  cat("\nAvailable Cell_annot labels:\n")
  print(sort(unique(na.omit(merged$Cell_annot))))
  stop(sprintf("NEUTROPHIL_LABEL '%s' not found in Cell_annot. Set it to one of the above.",
               NEUTROPHIL_LABEL))
}

# ── 2. SUBSET NEUTROPHILS ───────────────────────────────────
DefaultAssay(merged) <- "SCT"
Idents(merged) <- "Cell_annot"
neut <- subset(merged, idents = NEUTROPHIL_LABEL)

cat("\n=== Neutrophils per donor (Sample_Tag) ===\n")
print(table(neut$Sample_Tag, neut$Condition))

# ── 3. COHORT-LEVEL NEUTROPHIL DEG LIST (from script 07) ────
# Pull the neutrophil DEGs already computed cohort-wide so we test
# the SAME gene set for reproducibility.
all_degs <- read.csv("DEG_full_results_all_celltypes.csv",
                     stringsAsFactors = FALSE)

cohort_neut <- all_degs %>%
  filter(cell_type == NEUTROPHIL_LABEL) %>%
  mutate(cohort_dir = case_when(
    avg_log2FC >  LFC_CUT & p_val_adj < PADJ_CUT ~ "Up",
    avg_log2FC < -LFC_CUT & p_val_adj < PADJ_CUT ~ "Down",
    TRUE ~ "NS"
  )) %>%
  filter(cohort_dir != "NS")

cohort_deg_genes <- cohort_neut$gene
cat(sprintf("\nCohort-level neutrophil DEGs to test for consistency: %d\n",
            length(cohort_deg_genes)))

# ── 4. FOUR CROSS-DONOR PAIRINGS ────────────────────────────
# Each pairing = 1 diabetic donor vs 1 healthy donor.
# A robust DEG should recur, same direction, across pairings.
pairings <- expand.grid(diab = DIABETIC_TAGS,
                        healthy = HEALTHY_TAGS,
                        stringsAsFactors = FALSE)
pairings$label <- paste0(pairings$diab, "_vs_", pairings$healthy)

pairing_results <- list()

cat("\n=== Running per-pairing neutrophil DE (Diabetic donor vs Healthy donor) ===\n")
for (i in seq_len(nrow(pairings))) {
  d_tag <- pairings$diab[i]
  h_tag <- pairings$healthy[i]
  lab   <- pairings$label[i]
  
  cat(sprintf("  %-40s", lab))
  
  pair_cells <- subset(neut, subset = Sample_Tag %in% c(d_tag, h_tag))
  n_d <- sum(pair_cells$Sample_Tag == d_tag)
  n_h <- sum(pair_cells$Sample_Tag == h_tag)
  
  if (n_d < 10 || n_h < 10) {
    cat(sprintf("SKIPPED (diab=%d, healthy=%d)\n", n_d, n_h))
    next
  }
  
  # PrepSCTFindMarkers on the 2-donor subset so corrected counts
  # are recomputed for just this comparison
  pair_cells <- PrepSCTFindMarkers(pair_cells, verbose = FALSE)
  Idents(pair_cells) <- "Condition"
  
  tryCatch({
    mk <- FindMarkers(pair_cells,
                      ident.1  = "Diabetic",
                      ident.2  = "Healthy",
                      assay    = "SCT",
                      test.use = TEST_USE,
                      min.pct  = MIN_PCT,
                      verbose  = FALSE)
    mk$gene <- rownames(mk)
    mk$pairing <- lab
    mk$pair_dir <- case_when(
      mk$avg_log2FC >  LFC_CUT & mk$p_val_adj < PADJ_CUT ~ "Up",
      mk$avg_log2FC < -LFC_CUT & mk$p_val_adj < PADJ_CUT ~ "Down",
      TRUE ~ "NS"
    )
    pairing_results[[lab]] <- mk
    cat(sprintf("DEGs Up=%d Down=%d\n",
                sum(mk$pair_dir == "Up"), sum(mk$pair_dir == "Down")))
  }, error = function(e) cat(sprintf("ERROR: %s\n", e$message)))
}

pairing_all <- bind_rows(pairing_results)

# Export full per-pairing results (for per-pairing DEG-count figure in 08b)
write.csv(pairing_all, "neutrophil_DEG_per_pairing.csv", row.names = FALSE)

# Per-pairing DEG-count summary (Up / Down / total significant per pairing)
pairing_counts <- pairing_all %>%
  filter(pair_dir != "NS") %>%
  count(pairing, pair_dir) %>%
  tidyr::pivot_wider(names_from = pair_dir, values_from = n, values_fill = 0)
if (!"Up"   %in% names(pairing_counts)) pairing_counts$Up   <- 0
if (!"Down" %in% names(pairing_counts)) pairing_counts$Down <- 0
pairing_counts <- pairing_counts %>%
  mutate(Total = Up + Down) %>%
  arrange(desc(Total))
write.csv(pairing_counts, "neutrophil_DEG_per_pairing_counts.csv", row.names = FALSE)

# ── 5. CONSISTENCY SCORE PER COHORT DEG ─────────────────────
# For each cohort-level DEG, in how many of the 4 pairings is it
# recovered in the SAME direction as the cohort call?
n_pairings <- length(pairing_results)

consist <- lapply(cohort_deg_genes, function(g) {
  cohort_d <- cohort_neut$cohort_dir[cohort_neut$gene == g][1]
  sub <- pairing_all[pairing_all$gene == g, , drop = FALSE]
  # log2FC sign agreement (directional reproducibility, lenient)
  same_sign <- sum(sign(sub$avg_log2FC) ==
                     ifelse(cohort_d == "Up", 1, -1), na.rm = TRUE)
  # passes DEG threshold AND same direction (strict)
  same_dir_sig <- sum(sub$pair_dir == cohort_d, na.rm = TRUE)
  data.frame(
    gene              = g,
    cohort_direction  = cohort_d,
    n_pairings_tested = nrow(sub),
    n_same_sign       = same_sign,                       # /4 lenient
    n_same_dir_sig    = same_dir_sig,                    # /4 strict
    frac_same_sign    = round(same_sign / n_pairings, 2),
    stringsAsFactors  = FALSE
  )
}) %>% bind_rows() %>%
  arrange(desc(n_same_dir_sig), desc(n_same_sign))

# Classify reproducibility
consist <- consist %>%
  mutate(reproducibility = case_when(
    n_same_sign == n_pairings  ~ "Consistent (all pairings)",
    n_same_sign >= n_pairings - 1 ~ "Mostly consistent",
    n_same_sign <= 1           ~ "Outlier-driven / inconsistent",
    TRUE ~ "Partial"
  ))

cat("\n=== Reproducibility breakdown of cohort-level neutrophil DEGs ===\n")
print(table(consist$reproducibility))

write.csv(consist, "neutrophil_DEG_donor_consistency.csv", row.names = FALSE)

# ── 6. PER-DONOR MEAN EXPRESSION (so signal isn't one-donor) ─
# Mean SCT expression of each cohort DEG in each of the 4 donors.
sct_data <- tryCatch(
  GetAssayData(neut, assay = "SCT", layer = "data"),                  # Seurat v5
  error = function(e) GetAssayData(neut, assay = "SCT", slot = "data") # v4 fallback
)
genes_present <- intersect(cohort_deg_genes, rownames(sct_data))

donor_means <- sapply(c(DIABETIC_TAGS, HEALTHY_TAGS), function(tag) {
  cells <- colnames(neut)[neut$Sample_Tag == tag]
  Matrix::rowMeans(sct_data[genes_present, cells, drop = FALSE])
})
donor_means <- as.data.frame(donor_means) %>%
  rownames_to_column("gene")

write.csv(donor_means, "neutrophil_donor_mean_expression.csv", row.names = FALSE)

# ── 7. FIGURE A: per-donor expression heatmap of top DEGs ───
# z-scored across the 4 donors; donors annotated by condition.
top_genes <- consist %>%
  filter(reproducibility %in% c("Consistent (all pairings)", "Mostly consistent")) %>%
  slice_head(n = 50) %>% pull(gene)
top_genes <- intersect(top_genes, donor_means$gene)

if (length(top_genes) >= 2) {
  mat <- donor_means %>% filter(gene %in% top_genes) %>%
    column_to_rownames("gene") %>% as.matrix()
  mat <- mat[, c(HEALTHY_TAGS, DIABETIC_TAGS)]   # order columns H,H,D,D
  ann_col <- data.frame(
    Condition = c(rep("Healthy", length(HEALTHY_TAGS)),
                  rep("Diabetic", length(DIABETIC_TAGS))),
    row.names = colnames(mat)
  )
  hm <- pheatmap(mat,
                 scale = "row",
                 cluster_cols = FALSE,
                 cluster_rows = TRUE,
                 annotation_col = ann_col,
                 show_colnames = TRUE,
                 fontsize_row = 6,
                 main = "Per-donor expression of consistent neutrophil DEGs (z-scored)",
                 filename = "Fig_neutrophil_DEG_consistency_heatmap.pdf",
                 width = 5, height = max(6, length(top_genes) * 0.18))
  cat("\n✅ Saved: Fig_neutrophil_DEG_consistency_heatmap.pdf\n")
} else {
  cat("\n⚠ Not enough consistent genes to draw heatmap.\n")
}

# ── 8. FIGURE B: concordance of log2FC across pairings ──────
# Bar of how many cohort DEGs are recovered (same direction, sig)
# in 0,1,2,3,4 pairings — the headline reproducibility figure.
concord_tab <- consist %>%
  count(n_same_dir_sig) %>%
  mutate(n_same_dir_sig = factor(n_same_dir_sig,
                                 levels = 0:n_pairings))

p_concord <- ggplot(concord_tab,
                    aes(x = n_same_dir_sig, y = n)) +
  geom_col(fill = "#4575b4") +
  labs(x = sprintf("Number of donor pairings (of %d) recovering the DEG\n(same direction, adj P < 0.05)", n_pairings),
       y = "Number of cohort-level neutrophil DEGs",
       title = "Cross-donor reproducibility of the neutrophil DEG signature") +
  theme_classic(base_size = 11)

ggsave("Fig_neutrophil_donor_pairing_concordance.pdf",
       p_concord, width = 6, height = 4.2)
cat("✅ Saved: Fig_neutrophil_donor_pairing_concordance.pdf\n")

# ── 9. CONSOLE SUMMARY ──────────────────────────────────────
cat("\n==================================================\n")
cat("PER-DONOR CONSISTENCY SUMMARY\n")
cat("==================================================\n")
cat(sprintf("Cohort-level neutrophil DEGs tested : %d\n", nrow(consist)))
cat(sprintf("Pairings run                        : %d\n", n_pairings))
cat(sprintf("Consistent across ALL pairings      : %d (%.0f%%)\n",
            sum(consist$reproducibility == "Consistent (all pairings)"),
            100 * mean(consist$reproducibility == "Consistent (all pairings)")))
cat(sprintf("Outlier-driven / inconsistent       : %d (%.0f%%)\n",
            sum(consist$reproducibility == "Outlier-driven / inconsistent"),
            100 * mean(consist$reproducibility == "Outlier-driven / inconsistent")))
cat("\n✅ Section 07b complete\n")

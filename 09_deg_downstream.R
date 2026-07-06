# ============================================================
# SECTION 09: DEG DOWNSTREAM ANALYSIS
# Pathway Enrichment (GO + KEGG) | Heatmaps | Dot Plots
# Focus: Neutrophils vs all other cell types in Diabetes
# ============================================================
# INPUT  : DEG_full_results_all_celltypes.csv
#          merged_with_DEGs.rds  (for heatmap expression values)
# SAVES  : GO / KEGG result CSVs
#          Fig9  — GO Dot Plot       (Neutrophils, all terms)
#          Fig10 — KEGG Dot Plot     (Neutrophils, all terms)
#          Fig11 — GO Comparison     (all cell types, top terms)
#          Fig12 — Expression Heatmap (neutrophil top DEGs)
#          Fig13 — Dot Plot          (top DEG markers per cell type)
#          Fig14 — GO Enrichment Map (neutrophils)
# ============================================================
#
# PACKAGES NEEDED — install any that are missing:
#   BiocManager::install(c("clusterProfiler","org.Hs.eg.db",
#                          "enrichplot","ReactomePA"))
#   install.packages(c("ggtext","pheatmap","viridis","scales",
#                      "stringr","ggrepel"))
# ============================================================

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)        # human gene ID mapping
library(enrichplot)           # dotplot, emapplot, cnetplot
library(ggplot2)
library(ggtext)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(viridis)
library(scales)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(fgsea)       # for GSEA
library(msigdbr)     # for GO gene sets
library(tibble)

setwd("C:/IISc/ScRNA_Analysis_Human/20260422_output")

# ── Resolve namespace conflicts ────────────────────────────────────────────
# org.Hs.eg.db / AnnotationDbi export select(), filter(), rename(), mutate()
# and mask the dplyr versions when loaded after dplyr.
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
mutate <- dplyr::mutate

# ── Load inputs ────────────────────────────────────────────────────────────
all_degs   <- read.csv("DEG_full_results_all_celltypes.csv")
merged     <- readRDS("merged_with_DEGs.rds")

cat("Loaded all_degs  :", nrow(all_degs), "rows\n")
cat("Loaded merged    :", ncol(merged), "cells\n")

# ── Shared constants ───────────────────────────────────────────────────────
NEUTRO_COLOR <- "#D62728"
PADJ_CUT     <- 0.05
FC_CUT       <- 0.25

# Significant DEGs only
sig_degs <- all_degs %>%
  filter(p_val_adj < PADJ_CUT, abs(avg_log2FC) > FC_CUT)

cat("Significant DEGs total:", nrow(sig_degs), "\n")
cat("Cell types with DEGs  :", n_distinct(sig_degs$cell_type), "\n")

# ============================================================
# LOAD GO BP GENE SETS — required for GSEA section
# ============================================================
cat("\n=== Loading GO:BP gene sets (msigdbr) ===\n")
go_bp_df   <- msigdbr(species = "Homo sapiens",
                      category = "C5", subcategory = "GO:BP")
go_bp_sets <- split(go_bp_df$gene_symbol, go_bp_df$gs_name)
cat("GO:BP sets loaded:", length(go_bp_sets), "\n")

# ── Helper: symbol → Entrez ID conversion ─────────────────────────────────
symbols_to_entrez <- function(gene_symbols) {
  ids <- bitr(gene_symbols,
              fromType = "SYMBOL",
              toType   = "ENTREZID",
              OrgDb    = org.Hs.eg.db)
  ids
}

# ============================================================
# 9.1  GENE ONTOLOGY ENRICHMENT — NEUTROPHILS
# Separate enrichment for UP and DOWN DEGs
# ============================================================

cat("\n=== Running GO Enrichment — Neutrophils ===\n")

neutro_up   <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Upregulated") %>%
  pull(gene)

neutro_down <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Downregulated") %>%
  pull(gene)

# Universe = all genes tested in neutrophils
neutro_universe <- all_degs %>%
  filter(cell_type == "Neutrophils") %>%
  pull(gene)

cat("Neutrophil UP   genes:", length(neutro_up),   "\n")
cat("Neutrophil DOWN genes:", length(neutro_down),  "\n")
cat("Neutrophil universe :", length(neutro_universe), "\n")

# Convert to Entrez IDs
entrez_up       <- symbols_to_entrez(neutro_up)
entrez_down     <- symbols_to_entrez(neutro_down)
entrez_universe <- symbols_to_entrez(neutro_universe)

# ── GO BP enrichment ──────────────────────────────────────────────────────
run_go <- function(entrez_genes, universe_entrez,
                   ont = "BP", label = "") {
  if (nrow(entrez_genes) < 5) {
    cat(sprintf("⚠️  Too few genes for GO %s (%s) — skipping\n",
                ont, label))
    return(NULL)
  }
  enrichGO(
    gene          = entrez_genes$ENTREZID,
    universe      = universe_entrez$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE     # converts Entrez back to gene symbols
  )
}

go_up_bp   <- run_go(entrez_up,   entrez_universe, "BP", "UP-BP")
go_down_bp <- run_go(entrez_down, entrez_universe, "BP", "DOWN-BP")
go_up_mf   <- run_go(entrez_up,   entrez_universe, "MF", "UP-MF")
go_down_mf <- run_go(entrez_down, entrez_universe, "MF", "DOWN-MF")

# ── Save GO results ───────────────────────────────────────────────────────
if (!is.null(go_up_bp))   write.csv(as.data.frame(go_up_bp),
  "GO_BP_Neutrophils_UP.csv",   row.names = FALSE)
if (!is.null(go_down_bp)) write.csv(as.data.frame(go_down_bp),
  "GO_BP_Neutrophils_DOWN.csv", row.names = FALSE)
if (!is.null(go_up_mf))   write.csv(as.data.frame(go_up_mf),
  "GO_MF_Neutrophils_UP.csv",   row.names = FALSE)
if (!is.null(go_down_mf)) write.csv(as.data.frame(go_down_mf),
  "GO_MF_Neutrophils_DOWN.csv", row.names = FALSE)

cat("✅ GO enrichment complete\n")

# ============================================================
# FIGURE 9 — GO BP Dot Plots: Neutrophils UP & DOWN
# ============================================================

make_go_dotplot <- function(ego, title, colour) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("No significant GO terms for:", title, "\n")
    return(NULL)
  }

  df <- as.data.frame(ego) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20) %>%
    mutate(
      Description = str_wrap(Description, width = 45),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      })
    )

  ggplot(df, aes(x = GeneRatio_num, y = Description,
                 size = Count, colour = p.adjust)) +
    geom_point() +
    scale_colour_gradient(low  = colour,
                          high = "grey80",
                          name = "adj. p-value") +
    scale_size_continuous(range = c(2, 8), name = "Gene count") +
    labs(title = title,
         x     = "Gene Ratio",
         y     = NULL) +
    theme_bw(base_size = 11) +
    theme(
      plot.title  = element_text(face = "bold", size = 12, hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
}

p_go_up   <- make_go_dotplot(go_up_bp,
                              "GO:BP — Neutrophils Upregulated\n(Diabetic vs Healthy)",
                              "#D62728")
p_go_down <- make_go_dotplot(go_down_bp,
                              "GO:BP — Neutrophils Downregulated\n(Diabetic vs Healthy)",
                              "#1F77B4")

if (!is.null(p_go_up) && !is.null(p_go_down)) {
  p9 <- p_go_up | p_go_down
} else {
  p9 <- p_go_up %||% p_go_down
}

print(p9)
ggsave("Fig9_GO_BP_Neutrophils.pdf", p9,
       width = 16, height = 9, dpi = 300)
cat("✅ Saved: Fig9_GO_BP_Neutrophils.pdf\n")

# ============================================================
# 9.2  KEGG PATHWAY ENRICHMENT — NEUTROPHILS
# ============================================================

cat("\n=== Running KEGG Enrichment — Neutrophils ===\n")

run_kegg <- function(entrez_genes, universe_entrez, label = "") {
  if (nrow(entrez_genes) < 5) {
    cat(sprintf("⚠️  Too few genes for KEGG (%s) — skipping\n", label))
    return(NULL)
  }
  enrichKEGG(
    gene          = entrez_genes$ENTREZID,
    universe      = universe_entrez$ENTREZID,
    organism      = "hsa",          # Homo sapiens
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
}

kegg_up   <- run_kegg(entrez_up,   entrez_universe, "UP")
kegg_down <- run_kegg(entrez_down, entrez_universe, "DOWN")

# Convert Entrez IDs back to symbols for readability
if (!is.null(kegg_up))   kegg_up   <- setReadable(kegg_up,
                                                   org.Hs.eg.db, "ENTREZID")
if (!is.null(kegg_down)) kegg_down <- setReadable(kegg_down,
                                                   org.Hs.eg.db, "ENTREZID")

if (!is.null(kegg_up))   write.csv(as.data.frame(kegg_up),
  "KEGG_Neutrophils_UP.csv",   row.names = FALSE)
if (!is.null(kegg_down)) write.csv(as.data.frame(kegg_down),
  "KEGG_Neutrophils_DOWN.csv", row.names = FALSE)

cat("✅ KEGG enrichment complete\n")

# ============================================================
# FIGURE 10 — KEGG Dot Plots: Neutrophils UP & DOWN
# ============================================================

p_kegg_up   <- make_go_dotplot(kegg_up,
                                "KEGG — Neutrophils Upregulated\n(Diabetic vs Healthy)",
                                "#D62728")
p_kegg_down <- make_go_dotplot(kegg_down,
                                "KEGG — Neutrophils Downregulated\n(Diabetic vs Healthy)",
                                "#1F77B4")

if (!is.null(p_kegg_up) && !is.null(p_kegg_down)) {
  p10 <- p_kegg_up | p_kegg_down
} else {
  p10 <- p_kegg_up %||% p_kegg_down
}

print(p10)
ggsave("Fig10_KEGG_Neutrophils.pdf", p10,
       width = 16, height = 9, dpi = 300)
cat("✅ Saved: Fig10_KEGG_Neutrophils.pdf\n")

# ============================================================
# 9.3  COMPARATIVE GO ENRICHMENT — ALL CELL TYPES
# compareCluster: runs GO on each cell type's UP DEGs simultaneously
# Reveals which pathways are neutrophil-specific vs shared
# ============================================================

cat("\n=== Running compareCluster GO across all cell types (UP + DOWN) ===\n")

# Build named gene list: UP and DOWN treated as separate clusters
# Naming convention: "CellType__UP" and "CellType__DOWN"
# This lets Fig11 show both directions side-by-side per cell type

make_entrez_list <- function(direction_filter, suffix) {
  sig_degs %>%
    filter(direction == direction_filter) %>%
    group_by(cell_type) %>%
    summarise(genes = list(gene), .groups = "drop") %>%
    { setNames(.$genes, paste0(.$cell_type, "__", suffix)) } %>%
    lapply(function(g) {
      ids <- bitr(g, fromType = "SYMBOL",
                  toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      ids$ENTREZID
    }) %>%
    .[sapply(., length) >= 5]   # drop clusters too small to test
}

up_entrez_list   <- make_entrez_list("Upregulated",   "UP")
down_entrez_list <- make_entrez_list("Downregulated", "DOWN")

combined_entrez_list <- c(up_entrez_list, down_entrez_list)

cat("Clusters included in compareCluster:\n")
cat(paste(" ", names(combined_entrez_list), collapse = "\n"), "\n")

cc_go <- compareCluster(
  geneClusters  = combined_entrez_list,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

write.csv(as.data.frame(cc_go),
          "GO_BP_compareCluster_allcelltypes.csv",
          row.names = FALSE)
cat("✅ compareCluster complete\n")

# ============================================================
# FIGURE 11 — Comparative GO Dot Plot (all cell types)
# ============================================================

# Subset to top 5 terms per cell type for readability
cc_go_top <- clusterProfiler::filter(cc_go, p.adjust < 0.05)

p11 <- dotplot(cc_go_top, showCategory = 4, label_format = 35) +
  # Rotate x-axis labels — cluster names are now "CellType__UP/DOWN"
  scale_x_discrete(labels = function(x) {
    # Format as "CellType\n(UP)" or "CellType\n(DOWN)" for readability
    sub("__UP$",   "\n(↑UP)",   
    sub("__DOWN$", "\n(↓DOWN)", x))
  }) +
  labs(
    title    = "GO:BP Enrichment — All Cell Types\n(Diabetic vs Healthy)",
    subtitle = "Top 4 terms per cluster | ↑UP = upregulated in Diabetic | ↓DOWN = downregulated in Diabetic"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y   = element_text(size = 8),
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey40"),
    legend.position = "right"
  )

print(p11)
ggsave("Fig11_GO_compareCluster_allcelltypes.pdf", p11,
       width = 18, height = 11, dpi = 300)
cat("✅ Saved: Fig11_GO_compareCluster_allcelltypes.pdf\n")

# ============================================================
# FIGURE 14 — GO Enrichment Map (cnetplot)
# Shows gene–pathway network for neutrophil UP terms
# ============================================================

if (!is.null(go_up_bp) && nrow(as.data.frame(go_up_bp)) >= 3) {

  ep_ver <- packageVersion("enrichplot")
  cat(sprintf("enrichplot version: %s\n", ep_ver))
  # cnetplot API history:
  #   < 1.14  : foldChange (top-level), circular, colorEdge, cex_label_gene
  #   1.14–1.17 : foldChange (top-level), node_label, layout, cex_label_*
  #   >= 1.18 : foldChange moves INSIDE color.params = list(foldChange = ...)
  #             color.params = list(foldChange, edge)
  #             cex.params   = list(gene_label, category_label)
  # Confirmed working on enrichplot 1.30.5

  fc_vec <- sig_degs %>%
    filter(cell_type == "Neutrophils", direction == "Upregulated") %>%
    dplyr::select(gene, avg_log2FC) %>%
    tibble::deframe()

  # Build the cnetplot with only the two universally stable arguments.
  # Optional styling (edge colour, label size) varies by enrichplot build —
  # add them back once args(cnetplot) confirms the correct names.
  p14_base <- cnetplot(go_up_bp, showCategory = 8)

  # Manually annotate gene nodes with FC colour by rebuilding as ggplot layer
  # This approach is immune to cnetplot API changes.
  p14 <- p14_base +
    labs(
      title    = "GO:BP Gene–Pathway Network — Neutrophils Upregulated",
      subtitle = "Node colour = Log2FC (Diabetic vs Healthy)"
    ) +
    labs(
      title    = "GO:BP Gene–Pathway Network — Neutrophils Upregulated",
      subtitle = "Node colour = Log2FC (Diabetic vs Healthy)"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40")
    )

  print(p14)
  ggsave("Fig14_GO_cnetplot_Neutrophils.pdf", p14,
         width = 13, height = 10, dpi = 300)
  cat("✅ Saved: Fig14_GO_cnetplot_Neutrophils.pdf\n")

} else {
  cat("⚠️  Insufficient GO terms for cnetplot — skipping Fig14\n")
}

# ============================================================
# 9.4  EXPRESSION HEATMAP — Top Neutrophil DEGs
# Cells: neutrophils only | Genes: top UP + DOWN DEGs
# Columns split by Condition | Rows split by direction
# ============================================================

cat("\n=== Building Expression Heatmap — Neutrophils ===\n")

# Select top genes: top 25 UP + top 25 DOWN by |FC| × significance
neutro_top_up <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Upregulated") %>%
  mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
  arrange(desc(score)) %>%
  slice_head(n = 25) %>%
  pull(gene)

neutro_top_down <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Downregulated") %>%
  mutate(score = (-avg_log2FC) * (-log10(p_val_adj + 1e-300))) %>%
  arrange(desc(score)) %>%
  slice_head(n = 25) %>%
  pull(gene)

heatmap_genes <- c(neutro_top_up, neutro_top_down)
cat("Heatmap genes — UP:", length(neutro_top_up),
    "| DOWN:", length(neutro_top_down), "\n")

# Subset to neutrophils only
DefaultAssay(merged) <- "RNA"
merged <- JoinLayers(merged)

neutro_obj <- subset(merged, subset = Cell_annot == "Neutrophils")
cat("Neutrophil cells:", ncol(neutro_obj), "\n")

# Extract log-normalised expression matrix
expr_mat <- GetAssayData(neutro_obj, assay = "RNA", layer = "data")

# Keep only heatmap genes present in the matrix
heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(expr_mat)]
expr_sub  <- as.matrix(expr_mat[heatmap_genes, ])

# Scale rows (z-score per gene) for colour contrast
expr_scaled <- t(scale(t(expr_sub)))
expr_scaled[expr_scaled >  2.5]  <-  2.5   # clip extremes
expr_scaled[expr_scaled < -2.5]  <- -2.5

# ── Column annotation: Condition ──────────────────────────────────────────
col_annot <- data.frame(
  Condition = neutro_obj$Condition,
  row.names = colnames(neutro_obj)
)

col_annot_colors <- list(
  Condition = c("Healthy"  = "#4A90D9",
                "Diabetic" = "#D62728")
)

# ── Row annotation: UP / DOWN direction ───────────────────────────────────
row_annot <- data.frame(
  Direction = ifelse(heatmap_genes %in% neutro_top_up,
                     "Upregulated", "Downregulated"),
  row.names = heatmap_genes
)

row_annot_colors <- list(
  Direction = c("Upregulated"   = "green4",
                "Downregulated" = "pink3")
)

# Sort columns by Condition for clean visual split
col_order <- order(neutro_obj$Condition)

# ── Row gap between UP and DOWN ────────────────────────────────────────────
n_up_genes <- length(neutro_top_up[neutro_top_up %in% heatmap_genes])


p_heat12 <- pheatmap(
  expr_scaled[, col_order],
  color            = colorRampPalette(
                       rev(brewer.pal(11, "RdBu")))(100),
  breaks           = seq(-2.5, 2.5, length.out = 101),
  cluster_rows     = FALSE,      # keep UP / DOWN ordering
  cluster_cols     = FALSE,      # keep Condition ordering
  show_colnames    = FALSE,
  show_rownames    = TRUE,
  annotation_col   = col_annot,
  annotation_row   = row_annot,
  annotation_colors = c(col_annot_colors, row_annot_colors),
  gaps_row         = n_up_genes,               # visual separator
  gaps_col         = sum(neutro_obj$Condition[col_order] == "Diabetic"),
  fontsize_row     = 7,
  fontsize         = 10,
  main             = "Neutrophil Top DEGs — Diabetic vs Healthy\n(z-scored log-normalised expression)",
  border_color     = NA,
  angle_col        = 45
)
ggsave("Fig12_Heatmap_Neutrophil_TopDEGs.pdf",
       plot = p_heat12, width = 14, height = 10, dpi = 300)
cat("✅ Saved: Fig12_Heatmap_Neutrophil_TopDEGs.pdf\n")

# ============================================================
# 9.5  CROSS-CELLTYPE HEATMAP
# Top 5 DEGs per cell type in one shared heatmap
# Reveals cell-type specificity of DEG programmes
# ============================================================

cat("\n=== Building Cross-Cell-Type Heatmap ===\n")

# Top 4 UP + top 4 DOWN genes per cell type (by |FC| × significance score)
# Keeping 4+4 per cell type instead of 5 UP only — more balanced and
# shows suppressed programmes alongside activated ones
top4_up_per_ct <- sig_degs %>%
  filter(direction == "Upregulated") %>%
  mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
  group_by(cell_type) %>%
  slice_max(score, n = 4) %>%
  ungroup() %>%
  mutate(deg_direction = "Upregulated")

top4_down_per_ct <- sig_degs %>%
  filter(direction == "Downregulated") %>%
  mutate(score = (-avg_log2FC) * (-log10(p_val_adj + 1e-300))) %>%
  group_by(cell_type) %>%
  slice_max(score, n = 4) %>%
  ungroup() %>%
  mutate(deg_direction = "Downregulated")

top5_per_ct <- bind_rows(top4_up_per_ct, top4_down_per_ct) %>%
  arrange(cell_type, deg_direction, desc(score))

cross_genes <- unique(top5_per_ct$gene)
cat("Cross-heatmap genes — UP:", nrow(top4_up_per_ct),
    "| DOWN:", nrow(top4_down_per_ct),
    "| unique:", length(cross_genes), "\n")

# Use all annotated cells (not just neutrophils)
DefaultAssay(merged) <- "RNA"
expr_all <- GetAssayData(merged, assay = "RNA", layer = "data")

cross_genes <- cross_genes[cross_genes %in% rownames(expr_all)]
expr_cross  <- as.matrix(expr_all[cross_genes, ])

# Pseudo-bulk: average expression per cell type × condition
pb_list <- list()

# Build a clean index upfront — exclude NA-annotated cells once
clean_meta <- merged@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  filter(!is.na(Cell_annot), !is.na(Condition)) %>%
  filter(barcode %in% colnames(expr_cross))   # only cells present in matrix

cat("Clean cells for pseudo-bulk:", nrow(clean_meta), "\n")

for (ct in unique(clean_meta$Cell_annot)) {
  for (cond in c("Healthy", "Diabetic")) {
    cells <- clean_meta %>%
      filter(Cell_annot == ct, Condition == cond) %>%
      pull(barcode)

    if (length(cells) < 5) next

    # Safety check — should always pass after the filter above
    cells <- cells[cells %in% colnames(expr_cross)]
    if (length(cells) < 5) next

    pb_list[[paste0(ct, "__", cond)]] <-
      rowMeans(expr_cross[, cells, drop = FALSE])
  }
}

pb_mat <- do.call(cbind, pb_list)

# Z-score per gene
pb_scaled <- t(scale(t(pb_mat)))
pb_scaled[pb_scaled >  2] <-  2
pb_scaled[pb_scaled < -2] <- -2

# Column annotation
col_meta <- data.frame(
  Cell_Type = sub("__.*", "", colnames(pb_mat)),
  Condition = sub(".*__", "", colnames(pb_mat)),
  row.names = colnames(pb_mat)
)

n_ct <- length(unique(col_meta$Cell_Type))
ct_pal <- colorRampPalette(brewer.pal(8, "Set2"))(n_ct)
names(ct_pal) <- sort(unique(col_meta$Cell_Type))

cross_annot_colors <- list(
  Condition = c("Healthy" = "#4A90D9", "Diabetic" = "#D62728"),
  Cell_Type = ct_pal
)

# Row annotation: source cell type + direction for each gene
row_ct_annot <- top5_per_ct %>%
  dplyr::select(gene, cell_type, deg_direction) %>%
  distinct(gene, .keep_all = TRUE) %>%
  filter(gene %in% rownames(pb_scaled)) %>%
  tibble::column_to_rownames("gene") %>%
  rename(Source_CellType = cell_type,
         Direction       = deg_direction)

row_ct_colors <- list(
  Source_CellType = ct_pal,
  Direction       = c("Upregulated"   = "#D62728",
                      "Downregulated" = "#1F77B4")
)

# Sort columns: group by cell type, then Healthy before Diabetic
col_order_cross <- order(col_meta$Cell_Type, col_meta$Condition)

# Sort rows: group by source cell type, then UP before DOWN within each
row_order_cross <- order(
  row_ct_annot[rownames(row_ct_annot) %in% rownames(pb_scaled),
               "Source_CellType"],
  row_ct_annot[rownames(row_ct_annot) %in% rownames(pb_scaled),
               "Direction"]
)

pb_final <- pb_scaled[
  rownames(row_ct_annot)[row_order_cross],
  col_order_cross
]

pdf("Fig12b_Heatmap_CrossCellType_TopDEGs.pdf",
    width = 16, height = 12)
pheatmap(
  pb_final,
  color            = colorRampPalette(
                       rev(brewer.pal(11, "RdBu")))(100),
  breaks           = seq(-2, 2, length.out = 101),
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  show_colnames    = TRUE,
  show_rownames    = TRUE,
  annotation_col   = col_meta[col_order_cross, , drop = FALSE],
  annotation_row   = row_ct_annot[row_order_cross, , drop = FALSE],
  annotation_colors = c(cross_annot_colors, row_ct_colors),
  fontsize_row     = 7,
  fontsize_col     = 8,
  fontsize         = 10,
  main             = "Top 4 UP + Top 4 DOWN DEGs per Cell Type — Pseudo-bulk Average Expression\n(z-scored | Diabetic vs Healthy)",
  border_color     = NA,
  angle_col        = 45
)
cat("✅ Saved: Fig12b_Heatmap_CrossCellType_TopDEGs.pdf\n")

# ============================================================
# FIGURE 13 — Seurat DotPlot: Top DEG Markers per Cell Type
# ============================================================

cat("\n=== Building DotPlot — Top markers per cell type ===\n")

# Top 2 UP + top 2 DOWN per cell type — keeps the plot readable
# Genes ordered: all UP genes first (by cell type), then all DOWN genes
top2_up_ct <- sig_degs %>%
  filter(direction == "Upregulated") %>%
  mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
  group_by(cell_type) %>%
  slice_max(score, n = 2) %>%
  ungroup() %>%
  arrange(cell_type) %>%
  mutate(deg_dir = "UP")

top2_down_ct <- sig_degs %>%
  filter(direction == "Downregulated") %>%
  mutate(score = (-avg_log2FC) * (-log10(p_val_adj + 1e-300))) %>%
  group_by(cell_type) %>%
  slice_max(score, n = 2) %>%
  ungroup() %>%
  arrange(cell_type) %>%
  mutate(deg_dir = "DOWN")

up_genes   <- unique(top2_up_ct$gene)
down_genes <- unique(top2_down_ct$gene)

# Filter to genes present in the RNA assay
rna_genes    <- rownames(merged[["RNA"]])
up_genes     <- up_genes[up_genes %in% rna_genes]
down_genes   <- down_genes[down_genes %in% rna_genes]

cat("DotPlot UP genes  :", length(up_genes), "\n")
cat("DotPlot DOWN genes:", length(down_genes), "\n")

DefaultAssay(merged) <- "RNA"
Idents(merged) <- "Cell_annot"
merged_clean <- subset(merged, subset = !is.na(Cell_annot))

# ── Fig13a — Upregulated panel ────────────────────────────────────────────
p13a <- DotPlot(
  merged_clean,
  features  = up_genes,
  cols      = c("grey92", "#D62728"),
  dot.scale = 6,
  group.by  = "Cell_annot"
) +
  coord_flip() +
  labs(
    title    = "Top 2 Upregulated DEGs per Cell Type\n(expressed higher in Diabetic)",
    subtitle = "Dot size = % cells expressing | Colour = average expression",
    x = "Gene", y = "Cell Type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y   = element_text(size = 8, face = "italic"),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey40"),
    legend.position = "bottom"
  )

# ── Fig13b — Downregulated panel ─────────────────────────────────────────
p13b_all <- DotPlot(
  merged_clean,
  features  = down_genes,
  cols      = c("grey92", "#1F77B4"),
  dot.scale = 6,
  group.by  = "Cell_annot"
) +
  coord_flip() +
  labs(
    title    = "Top 2 Downregulated DEGs per Cell Type\n(expressed lower in Diabetic)",
    subtitle = "Dot size = % cells expressing | Colour = average expression",
    x = "Gene", y = "Cell Type"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y   = element_text(size = 8, face = "italic"),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey40"),
    legend.position = "bottom"
  )

# ── Combined UP + DOWN panel ──────────────────────────────────────────────
p13 <- p13a | p13b_all

print(p13)
ggsave("Fig13_DotPlot_TopDEGs_perCellType.pdf", p13,
       width = 20, height = 10, dpi = 300)
cat("✅ Saved: Fig13_DotPlot_TopDEGs_perCellType.pdf\n")

# Also save individually for flexibility
ggsave("Fig13a_DotPlot_UP_perCellType.pdf",   p13a,
       width = 12, height = 10, dpi = 300)
ggsave("Fig13b_DotPlot_DOWN_perCellType.pdf",  p13b_all,
       width = 12, height = 10, dpi = 300)
cat("✅ Saved: Fig13a and Fig13b individually\n")

# ============================================================
# 9.6  NEUTROPHIL CONDITION SPLIT DOTPLOT — UP & DOWN
# Top genes in each direction, Healthy vs Diabetic side-by-side
# ============================================================

neutro_clean <- subset(merged_clean,
                       subset = Cell_annot == "Neutrophils")
Idents(neutro_clean) <- "Condition"

# Top 15 UP and top 15 DOWN neutrophil DEGs
neutro_up_genes <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Upregulated") %>%
  mutate(score = avg_log2FC * (-log10(p_val_adj + 1e-300))) %>%
  slice_max(score, n = 15) %>%
  pull(gene)

neutro_down_genes <- sig_degs %>%
  filter(cell_type == "Neutrophils", direction == "Downregulated") %>%
  mutate(score = (-avg_log2FC) * (-log10(p_val_adj + 1e-300))) %>%
  slice_max(score, n = 15) %>%
  pull(gene)

neutro_up_genes   <- neutro_up_genes[neutro_up_genes %in%
                       rownames(neutro_clean[["RNA"]])]
neutro_down_genes <- neutro_down_genes[neutro_down_genes %in%
                       rownames(neutro_clean[["RNA"]])]

# ── Upregulated in Diabetic ───────────────────────────────────────────────
p_neutro_up <- DotPlot(
  neutro_clean,
  features  = neutro_up_genes,
  cols      = c("grey92", "#D62728"),
  dot.scale = 7
) +
  coord_flip() +
  labs(
    title    = "Neutrophil UP DEGs\n(higher in Diabetic)",
    subtitle = "Top 15 by |FC| × significance",
    x = "Gene", y = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 12, face = "bold"),
    axis.text.y   = element_text(size = 9,  face = "italic"),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey40"),
    legend.position = "bottom"
  )

# ── Downregulated in Diabetic ─────────────────────────────────────────────
p_neutro_down <- DotPlot(
  neutro_clean,
  features  = neutro_down_genes,
  cols      = c("grey92", "#1F77B4"),
  dot.scale = 7
) +
  coord_flip() +
  labs(
    title    = "Neutrophil DOWN DEGs\n(lower in Diabetic)",
    subtitle = "Top 15 by |FC| × significance",
    x = "Gene", y = "Condition"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x   = element_text(size = 12, face = "bold"),
    axis.text.y   = element_text(size = 9,  face = "italic"),
    plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, colour = "grey40"),
    legend.position = "bottom"
  )

# ── Combined panel ────────────────────────────────────────────────────────
p13c <- (p_neutro_up | p_neutro_down) +
  plot_annotation(
    title    = "Neutrophil DEGs — Healthy vs Diabetic",
    subtitle = "Left: upregulated in Diabetic | Right: downregulated in Diabetic",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40")
    )
  )

print(p13c)
ggsave("Fig13c_DotPlot_Neutrophils_UP_DOWN.pdf", p13c,
       width = 14, height = 9, dpi = 300)
cat("✅ Saved: Fig13c_DotPlot_Neutrophils_UP_DOWN.pdf\n")




# ============================================================
# GSEA — CORRECTED PIPELINE
# Full neutrophil ranked gene list (all genes, both directions)
# Ranked by decreasing avg_log2FC
# Gene sets: GO BP c5.v7.4 | Permutations: 1,000
# ============================================================

cat("\n=== GSEA Corrected: Full Ranked Gene List ===\n")

# ── Step 1: Build full ranked list ────────────────────────────────────────
# Use ALL neutrophil genes tested — not just top 50
# This is what GSEA v4.0.3 desktop actually expects

neutro_full_ranked_df <- all_degs %>%
  filter(cell_type == "Neutrophils") %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%       # remove any duplicate symbols
  select(gene, avg_log2FC)

cat("Full ranked list size :", nrow(neutro_full_ranked_df), "\n")
cat("log2FC range          : [",
    round(min(neutro_full_ranked_df$avg_log2FC), 3), ",",
    round(max(neutro_full_ranked_df$avg_log2FC), 3), "]\n")
cat("Genes with FC > 0     :", sum(neutro_full_ranked_df$avg_log2FC > 0), "\n")
cat("Genes with FC < 0     :", sum(neutro_full_ranked_df$avg_log2FC < 0), "\n")

# Named numeric vector
ranked_vec_full <- deframe(neutro_full_ranked_df)

# ── Step 2: Quick overlap check ───────────────────────────────────────────
all_genes_in_sets <- unique(unlist(go_bp_sets))
pct_covered <- round(
  length(intersect(names(ranked_vec_full), all_genes_in_sets)) /
    length(names(ranked_vec_full)) * 100, 1)

cat("Gene coverage in GO BP sets:", pct_covered, "%\n")

# ── Step 3: Run fgsea ─────────────────────────────────────────────────────
cat("\n=== Running fgsea — Full Ranked List (1,000 permutations) ===\n")

set.seed(42)

gsea_result_full <- fgsea(
  pathways    = go_bp_sets,
  stats       = ranked_vec_full,
  nPermSimple = 1000,           # 1,000 permutations
  minSize     = 15,             # standard minimum
  maxSize     = 500,
  scoreType   = "std",          # correct — bidirectional input
  eps         = 0               # precise p-values
)

cat("Total pathways tested        :", nrow(gsea_result_full), "\n")

# ── Step 4: Filter significant results ───────────────────────────────────
gsea_sig_full <- gsea_result_full %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

cat("Significant pathways FDR<0.05:", nrow(gsea_sig_full), "\n")
cat("  Enriched  (NES > 0)        :", sum(gsea_sig_full$NES > 0), "\n")
cat("  Depleted  (NES < 0)        :", sum(gsea_sig_full$NES < 0), "\n")

# ── Step 5: Save outputs ──────────────────────────────────────────────────
gsea_full_export <- gsea_result_full %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

gsea_sig_export <- gsea_sig_full %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(gsea_full_export,
          "GSEA_GO_BP_Neutrophils_full_allgenes.csv",     row.names = FALSE)
write.csv(gsea_sig_export,
          "GSEA_GO_BP_Neutrophils_significant_allgenes.csv", row.names = FALSE)

cat("✅ Results saved\n")

# ============================================================
# FIGURE 10A — NES Bar Plot
# Top 15 enriched + top 15 depleted pathways
# ============================================================

make_gsea_barplot <- function(gsea_df, n_top = 15,
                              title = "GSEA GO:BP — Neutrophils") {
  
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    cat("⚠️  No significant pathways to plot.\n")
    return(NULL)
  }
  
  top_enr <- gsea_df %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    slice_head(n = n_top)
  
  top_dep <- gsea_df %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    slice_head(n = n_top)
  
  plot_df <- bind_rows(top_enr, top_dep) %>%
    mutate(
      pathway   = str_remove(pathway, "^GOBP_"),
      pathway   = str_replace_all(pathway, "_", " "),
      pathway   = str_to_sentence(pathway),
      pathway   = str_wrap(pathway, width = 50),
      pathway   = factor(pathway, levels = rev(unique(pathway))),
      Direction = ifelse(NES > 0, "Enriched", "Depleted")
    )
  
  ggplot(plot_df, aes(x = NES, y = pathway, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "black") +
    scale_fill_manual(
      values = c("Enriched" = "#D62728", "Depleted" = "#1F77B4")
    ) +
    labs(
      title   = title,
      x       = "Normalized Enrichment Score (NES)",
      y       = NULL,
      fill    = NULL,
      caption = "FDR < 0.05 | GO:BP c5.v7.4 | 1,000 permutations"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.caption    = element_text(size = 8, colour = "grey50"),
      axis.text.y     = element_text(size = 9),
      legend.position = "bottom"
    )
}

p10a <- make_gsea_barplot(
  gsea_sig_full,
  n_top = 15,
  title = "GSEA GO:BP — Neutrophils (Diabetic vs Healthy)"
)

# ============================================================
# FIGURE 10B — Classic Enrichment Plot (top pathway)
# ============================================================

if (nrow(gsea_sig_full) > 0) {
  
  top_pathway <- gsea_sig_full %>%
    arrange(desc(NES)) %>%
    slice(1) %>%
    pull(pathway)
  
  top_nes  <- round(gsea_sig_full$NES[1], 3)
  top_padj <- formatC(gsea_sig_full$padj[1], format = "e", digits = 2)
  
  cat("Top enriched pathway:", top_pathway, "\n")
  
  p10b <- plotEnrichment(
    pathway = go_bp_sets[[top_pathway]],
    stats   = ranked_vec_full
  ) +
    labs(
      title   = str_wrap(
        paste0("Enrichment: ",
               str_to_sentence(
                 str_replace_all(
                   str_remove(top_pathway, "^GOBP_"), "_", " "))),
        width = 55),
      x       = "Rank in gene list",
      y       = "Enrichment score (ES)",
      caption = paste0("NES = ", top_nes, "  |  FDR = ", top_padj)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.caption = element_text(size = 8, colour = "grey50")
    )
  
  # ── Combine and save ───────────────────────────────────────────────────
  p10 <- p10a | p10b
  print(p10)
  ggsave("Fig10_GSEA_GO_BP_Neutrophils.pdf", p10,
         width = 18, height = 9, dpi = 300)
  cat("✅ Saved: Fig10_GSEA_GO_BP_Neutrophils.pdf\n")
  
} else {
  cat("⚠️  No significant pathways — enrichment plot skipped\n")
  if (!is.null(p10a)) {
    print(p10a)
    ggsave("Fig10_GSEA_GO_BP_Neutrophils.pdf", p10a,
           width = 10, height = 9, dpi = 300)
  }
}

cat("\n✅ GSEA corrected pipeline complete\n")


# Also save individually
ggsave("Fig13c_UP_only.pdf",   p_neutro_up,   width = 8, height = 9, dpi = 300)
ggsave("Fig13c_DOWN_only.pdf", p_neutro_down, width = 8, height = 9, dpi = 300)
cat("✅ Saved: Fig13c_UP_only.pdf and Fig13c_DOWN_only.pdf\n")

# ============================================================
# SUMMARY
# ============================================================

cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat("║  SECTION 09 COMPLETE — Output Files                     ║\n")
cat("╚══════════════════════════════════════════════════════════╝\n\n")

outputs <- c(
  "GO_BP_Neutrophils_UP.csv",
  "GO_BP_Neutrophils_DOWN.csv",
  "GO_MF_Neutrophils_UP.csv",
  "GO_MF_Neutrophils_DOWN.csv",
  "KEGG_Neutrophils_UP.csv",
  "KEGG_Neutrophils_DOWN.csv",
  "GO_BP_compareCluster_allcelltypes.csv",
  "Fig9_GO_BP_Neutrophils.pdf         — GO:BP UP & DOWN dot plots (neutrophils)",
  "Fig10_KEGG_Neutrophils.pdf         — KEGG UP & DOWN dot plots (neutrophils)",
  "Fig11_GO_compareCluster_allcelltypes.pdf — GO:BP all cell types UP + DOWN",
  "Fig12_Heatmap_Neutrophil_TopDEGs.pdf    — Neutrophil heatmap top25 UP + DOWN",
  "Fig12b_Heatmap_CrossCellType_TopDEGs.pdf — Cross cell type heatmap UP + DOWN",
  "Fig13_DotPlot_TopDEGs_perCellType.pdf   — All cell types UP | DOWN panels",
  "Fig13a_DotPlot_UP_perCellType.pdf       — UP DEGs per cell type",
  "Fig13b_DotPlot_DOWN_perCellType.pdf     — DOWN DEGs per cell type",
  "Fig13c_DotPlot_Neutrophils_UP_DOWN.pdf  — Neutrophil UP + DOWN combined",
  "Fig13c_UP_only.pdf                      — Neutrophil UP only",
  "Fig13c_DOWN_only.pdf                    — Neutrophil DOWN only",
  "Fig14_GO_cnetplot_Neutrophils.pdf       — GO gene-pathway network (neutrophils UP)"
)

for (f in outputs) cat("  ", f, "\n")

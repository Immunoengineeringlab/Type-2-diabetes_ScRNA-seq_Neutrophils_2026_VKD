# ============================================================
# SECTION 08: DEG VISUALISATION — Neutrophils Most Altered
# ============================================================
# INPUT  : DEG_summary_all_celltypes.csv
#          DEG_full_results_all_celltypes.csv
# OUTPUT : 8 publication-quality figures
# SAVES  : Fig1–Fig8 PDFs
#
# NOTE: requires ggtext — install if needed:
#   install.packages("ggtext")
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(scales)
library(ggtext)

# ── Load from CSVs (no Seurat object needed for this section) ─────────────
summary_df <- read.csv("DEG_summary_all_celltypes.csv")
all_degs   <- read.csv("DEG_full_results_all_celltypes.csv")

cat("Loaded summary_df:", nrow(summary_df), "cell types\n")
cat("Loaded all_degs  :", nrow(all_degs),   "gene rows\n")

# ============================================================
# SHARED SETUP
# ============================================================

NEUTRO_COLOR <- "#D62728"

cell_types_ordered <- summary_df %>%
  arrange(desc(Total_DEGs)) %>%
  pull(Cell_Type)

n_types       <- length(cell_types_ordered)
base_pal      <- colorRampPalette(brewer.pal(8, "Set2"))(n_types)
names(base_pal) <- cell_types_ordered

highlight_pal <- ifelse(names(base_pal) == "Neutrophils",
                        NEUTRO_COLOR, "#BBBBBB")
names(highlight_pal) <- names(base_pal)

# Factor: most DEGs at top
summary_df$Cell_Type <- factor(summary_df$Cell_Type,
                               levels = rev(cell_types_ordered))

# Helper: HTML-styled axis labels (bold red for neutrophils)
make_html_labels <- function(lvls,
                             highlight = "Neutrophils",
                             hi_col    = NEUTRO_COLOR,
                             other_col = "grey25") {
  ifelse(lvls == highlight,
         sprintf("<b style='color:%s'>%s</b>", hi_col, lvls),
         sprintf("<span style='color:%s'>%s</span>", other_col, lvls))
}

# ============================================================
# FIGURE 1 — Lollipop: Total DEGs per Cell Type
# ============================================================

lollipop_labels <- make_html_labels(levels(summary_df$Cell_Type))

p1 <- ggplot(summary_df,
             aes(x = Cell_Type, y = Total_DEGs, colour = Cell_Type)) +
  geom_segment(aes(xend = Cell_Type, y = 0, yend = Total_DEGs),
               colour = "grey60", linewidth = 0.8) +
  geom_point(size = 5, shape = 21, fill = "white", stroke = 1.5) +
  geom_point(size = 3.5) +
  geom_text(aes(label = Total_DEGs),
            hjust = -0.55, size = 3.8, fontface = "bold",
            colour = "grey20") +
  scale_colour_manual(values = highlight_pal) +
  scale_x_discrete(labels = lollipop_labels) +
  coord_flip(clip = "off") +
  expand_limits(y = max(summary_df$Total_DEGs) * 1.15) +
  labs(title    = "Total Significant DEGs per Cell Type\n(Diabetic vs Healthy)",
       subtitle = "Adjusted p < 0.05 & |log2FC| > 0.25",
       x = NULL, y = "Number of DEGs") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.text.y        = element_markdown(size = 11),
    panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed")
  )

print(p1)
ggsave("Fig1_DEG_lollipop_per_celltype.pdf", p1,
       width = 8, height = 5.5, dpi = 300)
cat("✅ Saved: Fig1_DEG_lollipop_per_celltype.pdf\n")

# ============================================================
# FIGURE 2 — Stacked Bar: UP vs DOWN per Cell Type
# ============================================================

bar_df <- summary_df %>%
  select(Cell_Type, UP, DOWN) %>%
  pivot_longer(c(UP, DOWN), names_to = "Direction", values_to = "Count") %>%
  mutate(Direction    = factor(Direction, levels = c("UP", "DOWN")),
         Count_signed = ifelse(Direction == "DOWN", -Count, Count))

bar_labels  <- make_html_labels(levels(summary_df$Cell_Type))
neutro_pos  <- which(levels(summary_df$Cell_Type) == "Neutrophils")

p2 <- ggplot(bar_df,
             aes(x = Cell_Type, y = Count_signed, fill = Direction)) +
  annotate("rect",
           xmin = neutro_pos - 0.45, xmax = neutro_pos + 0.45,
           ymin = -Inf, ymax = Inf,
           fill = NEUTRO_COLOR, alpha = 0.07) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.5) +
  scale_fill_manual(
    values = c("UP" = "#E05A5A", "DOWN" = "#4A90D9"),
    labels = c("UP" = "↑ Upregulated in Diabetic",
               "DOWN" = "↓ Downregulated in Diabetic")
  ) +
  scale_y_continuous(labels = function(x) abs(x),
                     breaks = pretty_breaks(n = 6)) +
  scale_x_discrete(labels = bar_labels) +
  coord_flip() +
  labs(title    = "Directional DEGs per Cell Type\n(Diabetic vs Healthy)",
       subtitle = "Bar length = number of DEGs in each direction",
       x = NULL, y = "Number of DEGs", fill = NULL) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 11),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.text.y     = element_markdown(size = 11)
  )

print(p2)
ggsave("Fig2_DEG_directional_bar.pdf", p2,
       width = 8, height = 5.5, dpi = 300)
cat("✅ Saved: Fig2_DEG_directional_bar.pdf\n")

# ============================================================
# FIGURE 3 — Bubble: UP vs DOWN 2D landscape
# ============================================================

p3 <- ggplot(summary_df,
             aes(x = UP, y = DOWN, size = Total_DEGs,
                 colour = Cell_Type, label = Cell_Type)) +
  geom_point(alpha = 0.85) +
  geom_text_repel(
    aes(fontface = ifelse(Cell_Type == "Neutrophils", "bold", "plain")),
    size          = 3.8,
    colour        = ifelse(summary_df$Cell_Type == "Neutrophils",
                           NEUTRO_COLOR, "grey30"),
    box.padding   = 0.5,
    point.padding = 0.3,
    max.overlaps  = 20
  ) +
  scale_colour_manual(values = highlight_pal) +
  scale_size_continuous(range = c(3, 14), name = "Total DEGs") +
  labs(title    = "DEG Landscape Across Cell Types",
       subtitle = "Bubble size = Total DEGs | Axes = directional breakdown",
       x = "↑ Upregulated in Diabetic",
       y = "↓ Downregulated in Diabetic") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, colour = "grey40")
  ) +
  guides(colour = "none")

print(p3)
ggsave("Fig3_DEG_bubble_UP_vs_DOWN.pdf", p3,
       width = 8, height = 6.5, dpi = 300)
cat("✅ Saved: Fig3_DEG_bubble_UP_vs_DOWN.pdf\n")

# ============================================================
# FIGURE 4 — Volcano: Neutrophils Only
# ============================================================

neutro_degs <- all_degs %>%
  filter(cell_type == "Neutrophils") %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    direction = dplyr::case_when(
      avg_log2FC >  0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

top_up   <- neutro_degs %>% filter(direction == "Up") %>%
  slice_max(avg_log2FC * neg_log10_padj, n = 15)
top_down <- neutro_degs %>% filter(direction == "Down") %>%
  slice_max(-avg_log2FC * neg_log10_padj, n = 15)
label_genes <- bind_rows(top_up, top_down)

p4 <- ggplot(neutro_degs,
             aes(x = avg_log2FC, y = neg_log10_padj,
                 colour = direction)) +
  geom_point(data   = filter(neutro_degs, direction == "NS"),
             size = 1.2, alpha = 0.35, colour = "grey70") +
  geom_point(data   = filter(neutro_degs, direction != "NS"),
             size = 1.8, alpha = 0.80) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_text_repel(data          = label_genes,
                  aes(label     = gene),
                  size          = 3.2,
                  fontface      = "italic",
                  max.overlaps  = 30,
                  box.padding   = 0.4,
                  segment.colour = "grey50",
                  segment.size  = 0.3) +
  scale_colour_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey70"),
    labels = c("Up" = "↑ Upregulated", "Down" = "↓ Downregulated",
               "NS" = "Not Significant"),
    name   = NULL
  ) +
  labs(
    title    = "Volcano Plot — Neutrophils\n(Diabetic vs Healthy)",
    subtitle = paste0(sum(neutro_degs$direction == "Up"),
                      " upregulated | ",
                      sum(neutro_degs$direction == "Down"),
                      " downregulated"),
    x = "Average Log2 Fold Change",
    y = "-log10(Adjusted p-value)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, colour = "grey40")
  )

print(p4)
ggsave("Fig4_Volcano_Neutrophils.pdf", p4,
       width = 8, height = 7, dpi = 300)
cat("✅ Saved: Fig4_Volcano_Neutrophils.pdf\n")


# ============================================================
# FIGURE 4 — Volcano: Neutrophils, labelled by
# PER-DONOR-REPRODUCIBLE genes (not top-by-fold-change).
# ============================================================
# ---- Curated reproducible neutrophil genes (from 07b consistency) ----
# These are the genes shown in the per-donor consistency heatmap
REPRODUCIBLE_GENES <- c(
  # Up in T2D (reproducible)
  "PADI2","CXCL8","EGR3","TAGAP","GPR65","TRIB1","SAT1","ANP32A",
  # Down in T2D (reproducible)
  "CLEC12A","IRS2","SOCS3","RNASE6","CXCR1","CXCR2","CEBPB","CEBPD",
  "LYZ","THBD","TNFRSF1B","PILRA","VNN2","ASAH1","CHI3L1","S100P",
  "PECAM1","FGL2","CST3"
)

neutro_degs <- all_degs %>%
  filter(cell_type == "Neutrophils") %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    direction = dplyr::case_when(
      avg_log2FC >  0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    reproducible = gene %in% REPRODUCIBLE_GENES
  )

# Label only reproducible genes that are also significant
label_genes <- neutro_degs %>%
  filter(reproducible, direction != "NS")

cat(sprintf("Labelling %d reproducible significant genes on the volcano\n",
            nrow(label_genes)))
missing_repro <- setdiff(REPRODUCIBLE_GENES, neutro_degs$gene[neutro_degs$direction != "NS"])
if (length(missing_repro) > 0)
  cat("Reproducible genes not significant at cohort level (not labelled):\n  ",
      paste(missing_repro, collapse = ", "), "\n")

p4 <- ggplot(neutro_degs,
             aes(x = avg_log2FC, y = neg_log10_padj, colour = direction)) +
  # all non-significant points, faint
  geom_point(data = filter(neutro_degs, direction == "NS"),
             size = 1.2, alpha = 0.30, colour = "grey75") +
  # significant but not-reproducible points: plotted, coloured, slightly faded
  geom_point(data = filter(neutro_degs, direction != "NS", !reproducible),
             size = 1.6, alpha = 0.55) +
  # reproducible significant points: emphasised (outlined, full opacity)
  geom_point(data = filter(neutro_degs, direction != "NS", reproducible),
             size = 2.6, alpha = 0.95, shape = 21, stroke = 0.6,
             aes(fill = direction), colour = "grey20") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_text_repel(data = label_genes,
                  aes(label = gene),
                  size = 3.4, fontface = "italic",
                  max.overlaps = 40, box.padding = 0.5,
                  segment.colour = "grey50", segment.size = 0.3,
                  min.segment.length = 0) +
  scale_colour_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey75"),
    labels = c("Up" = "↑ Up in T2D", "Down" = "↓ Down in T2D",
               "NS" = "Not significant"),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4"), guide = "none"
  ) +
  labs(
    title    = "Volcano Plot — Neutrophils (Diabetic vs Healthy)",
    subtitle = paste0(sum(neutro_degs$direction == "Up"), " up | ",
                      sum(neutro_degs$direction == "Down"),
                      " down  •  labelled = reproducible across donors"),
    x = "Average log2 fold change",
    y = "-log10(adjusted P)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 9.5, hjust = 0.5, colour = "grey40")
  )

print(p4)
ggsave("Fig4_Volcano_Neutrophils.pdf", p4,
       width = 8, height = 7, dpi = 300)
cat("✅ Saved: Fig4_Volcano_Neutrophils.pdf (reproducible-gene labelling)\n")
# ============================================================
# FIGURE 5 — Multi-volcano grid (all cell types)
# ============================================================

all_degs_plot <- all_degs %>%
  mutate(
    neg_log10_padj = -log10(p_val_adj + 1e-300),
    direction = dplyr::case_when(
      avg_log2FC >  0.25 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    cell_type_label = ifelse(
      cell_type == "Neutrophils",
      sprintf("<b style='color:%s'>%s</b>", NEUTRO_COLOR, cell_type),
      cell_type
    )
  )

p5 <- ggplot(all_degs_plot,
             aes(x = avg_log2FC, y = neg_log10_padj,
                 colour = direction)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.25, 0.25),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  facet_wrap(~ cell_type_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey80"),
    name   = NULL
  ) +
  labs(title    = "Volcano Plots — All Cell Types\n(Diabetic vs Healthy)",
       subtitle = "Dashed lines: padj = 0.05 & |log2FC| = 0.25",
       x = "Log2 Fold Change", y = "-log₁₀(adj. p)") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_markdown(size = 10),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, colour = "grey40"),
    panel.grid      = element_blank()
  )

print(p5)
ggsave("Fig5_Volcano_AllCellTypes.pdf", p5,
       width = 13, height = 9, dpi = 300)
cat("✅ Saved: Fig5_Volcano_AllCellTypes.pdf\n")

# ============================================================
# FIGURE 6 — % DEGs normalised by genes tested
# ============================================================

genes_tested <- all_degs %>%
  group_by(cell_type) %>%
  summarise(genes_tested = n(), .groups = "drop")

summary_norm <- summary_df %>%
  left_join(genes_tested, by = c("Cell_Type" = "cell_type")) %>%
  mutate(pct_DEGs = Total_DEGs / genes_tested * 100) %>%
  arrange(desc(pct_DEGs)) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = rev(as.character(Cell_Type))))

norm_labels <- make_html_labels(levels(summary_norm$Cell_Type))

p6 <- ggplot(summary_norm,
             aes(x = Cell_Type, y = pct_DEGs, fill = Cell_Type)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(pct_DEGs, 1), "%")),
            hjust = -0.15, size = 3.5, fontface = "bold",
            colour = "grey20") +
  scale_fill_manual(values = highlight_pal) +
  scale_x_discrete(labels = norm_labels) +
  coord_flip(clip = "off") +
  expand_limits(y = max(summary_norm$pct_DEGs, na.rm = TRUE) * 1.2) +
  labs(title    = "DEGs as % of Tested Genes per Cell Type",
       subtitle = "Corrects for transcriptome depth differences",
       x = NULL, y = "% Significant DEGs") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.text.y        = element_markdown(size = 11),
    panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed")
  )

print(p6)
ggsave("Fig6_DEG_pct_normalised.pdf", p6,
       width = 8, height = 5.5, dpi = 300)
cat("✅ Saved: Fig6_DEG_pct_normalised.pdf\n")

# ============================================================
# FIGURE 7 — Combined panel (Fig 1 + Fig 4)
# ============================================================

combined <- (p1 / p4) +
  plot_annotation(
    title    = "Neutrophils Are the Most Transcriptionally Altered\nCell Type in Diabetes",
    subtitle = "Peripheral blood WBCs | BD Rhapsody | SCT DEG (Wilcoxon)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey40")
    )
  ) +
  plot_layout(heights = c(1, 1.2))

print(combined)
ggsave("Fig7_Combined_Panel.pdf", combined,
       width = 9, height = 13, dpi = 300)
cat("✅ Saved: Fig7_Combined_Panel.pdf\n")

# ============================================================
# FIGURE 8 — |log2FC| Distribution per Cell Type (Boxplot)
# Shows neutrophils have larger magnitude changes, not just more DEGs
# ============================================================

sig_degs <- all_degs %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
  mutate(cell_type = factor(cell_type, levels = cell_types_ordered))

fc_labels <- make_html_labels(rev(cell_types_ordered))

p8 <- ggplot(sig_degs,
             aes(x = cell_type, y = abs(avg_log2FC), fill = cell_type)) +
  geom_boxplot(outlier.size  = 0.8, outlier.alpha = 0.5,
               width = 0.65, colour = "grey30") +
  scale_fill_manual(values = highlight_pal) +
  scale_x_discrete(limits = rev(cell_types_ordered),
                   labels = fc_labels) +
  coord_flip() +
  labs(title    = "|Log2 FC| Distribution of Significant DEGs",
       subtitle = "Shows magnitude of change, not just count",
       x = NULL, y = "|Average Log2 Fold Change|") +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    axis.text.y        = element_markdown(size = 11),
    panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed")
  )

print(p8)
ggsave("Fig8_FC_distribution_per_celltype.pdf", p8,
       width = 8, height = 5.5, dpi = 300)
cat("✅ Saved: Fig8_FC_distribution_per_celltype.pdf\n")

# ============================================================
# SUMMARY
# ============================================================

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║   DEG SUMMARY — sorted by Total DEGs                ║\n")
cat("╚══════════════════════════════════════════════════════╝\n\n")

summary_df %>%
  arrange(desc(Total_DEGs)) %>%
  mutate(`%UP`   = round(UP   / Total_DEGs * 100, 1),
         `%DOWN` = round(DOWN / Total_DEGs * 100, 1)) %>%
  as.data.frame() %>%
  print()

cat("\n✅ Section 08 complete — all figures saved\n")
cat("Output files:\n")
cat("  Fig1_DEG_lollipop_per_celltype.pdf\n")
cat("  Fig2_DEG_directional_bar.pdf\n")
cat("  Fig3_DEG_bubble_UP_vs_DOWN.pdf\n")
cat("  Fig4_Volcano_Neutrophils.pdf\n")
cat("  Fig5_Volcano_AllCellTypes.pdf\n")
cat("  Fig6_DEG_pct_normalised.pdf\n")
cat("  Fig7_Combined_Panel.pdf\n")
cat("  Fig8_FC_distribution_per_celltype.pdf\n")
# ============================================================
# SECTION 08b: PER-DONOR DEG VISUALISATION (addendum to Section 08)
# ============================================================
# Renders the per-donor reproducibility figures in the same house
# style as Section 08. Reads the CSV outputs of Section 07b
# (no Seurat object needed).
#
# INPUT  : neutrophil_DEG_donor_consistency.csv      (from 07b)
#          neutrophil_donor_mean_expression.csv       (from 07b)
# OUTPUT : Fig9_donor_pairing_concordance.pdf
#          Fig10_donor_expression_heatmap.pdf
#          Fig11_donor_topDEG_dotplot.pdf
#
# NOTE: run AFTER Section 07b (which writes the two CSVs above) and
#       AFTER Section 08 (which defines NEUTRO_COLOR etc.). If running
#       08b standalone, the SHARED SETUP block below re-defines them.
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggtext)
library(pheatmap)
library(scales)

# ── SHARED SETUP (re-defined so 08b can run standalone) ─────
if (!exists("NEUTRO_COLOR")) NEUTRO_COLOR <- "#D62728"
HEALTHY_TAGS  <- c("SampleTag04_hs", "SampleTag12_hs")
DIABETIC_TAGS <- c("SampleTag01_hs", "SampleTag05_hs")

# ── LOAD 07b OUTPUTS ────────────────────────────────────────
stopifnot(file.exists("neutrophil_DEG_donor_consistency.csv"))
stopifnot(file.exists("neutrophil_donor_mean_expression.csv"))

consist     <- read.csv("neutrophil_DEG_donor_consistency.csv",
                        stringsAsFactors = FALSE)
donor_means <- read.csv("neutrophil_donor_mean_expression.csv",
                        stringsAsFactors = FALSE)

cat("Loaded consistency table:", nrow(consist), "DEGs\n")
cat("Loaded donor means      :", nrow(donor_means), "DEGs x",
    ncol(donor_means) - 1, "donors\n")

n_pairings <- max(consist$n_pairings_tested, na.rm = TRUE)

# ============================================================
# FIGURE 9 — Cross-donor concordance bar (house style)
# How many of the N pairings recover each cohort DEG (strict:
# same direction AND adj p < 0.05)
# ============================================================

concord_tab <- consist %>%
  count(n_same_dir_sig) %>%
  tidyr::complete(n_same_dir_sig = 0:n_pairings, fill = list(n = 0)) %>%
  mutate(n_same_dir_sig = factor(n_same_dir_sig, levels = 0:n_pairings))

p9 <- ggplot(concord_tab, aes(x = n_same_dir_sig, y = n)) +
  geom_col(width = 0.7, fill = NEUTRO_COLOR, alpha = 0.85) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.8,
            fontface = "bold", colour = "grey20") +
  expand_limits(y = max(concord_tab$n) * 1.12) +
  labs(
    title    = "Cross-Donor Reproducibility of the Neutrophil DEG Signature",
    subtitle = sprintf("Strict: same direction & adj p < 0.05 in each of %d donor pairings", n_pairings),
    x = sprintf("Number of donor pairings (of %d) recovering the DEG", n_pairings),
    y = "Number of cohort-level neutrophil DEGs"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed")
  )

print(p9)
ggsave("Fig9_donor_pairing_concordance.pdf", p9,
       width = 8, height = 5.5, dpi = 300)
cat("\u2705 Saved: Fig9_donor_pairing_concordance.pdf\n")

# ============================================================
# FIGURE 10 — Per-donor expression heatmap of top consistent DEGs
# z-scored across the 4 donors; columns ordered Healthy, Diabetic
# ============================================================

# Pick top consistent DEGs: prefer the lenient sign-agreement metric
# (robust to per-pairing power loss), then break ties by strict count.
top_genes <- consist %>%
  arrange(desc(n_same_sign), desc(n_same_dir_sig)) %>%
  filter(n_same_sign >= max(1, n_pairings - 1)) %>%   # consistent direction
  slice_head(n = 50) %>%
  pull(gene)
top_genes <- intersect(top_genes, donor_means$gene)

donor_cols <- intersect(c(HEALTHY_TAGS, DIABETIC_TAGS), colnames(donor_means))

if (length(top_genes) >= 2 && length(donor_cols) >= 2) {
  mat <- donor_means %>%
    filter(gene %in% top_genes) %>%
    column_to_rownames("gene") %>%
    select(all_of(donor_cols)) %>%
    as.matrix()
  
  ann_col <- data.frame(
    Condition = ifelse(colnames(mat) %in% DIABETIC_TAGS, "Diabetic", "Healthy"),
    row.names = colnames(mat)
  )
  ann_colors <- list(Condition = c(Healthy = "#9a9aa0", Diabetic = NEUTRO_COLOR))
  
  # pheatmap writes its own PDF via filename= (no dev.off needed)
  pheatmap(
    mat,
    scale            = "row",
    cluster_cols     = FALSE,
    cluster_rows     = TRUE,
    annotation_col   = ann_col,
    annotation_colors = ann_colors,
    show_colnames    = TRUE,
    fontsize_row     = 6,
    color            = colorRampPalette(c("#1F77B4", "white", NEUTRO_COLOR))(100),
    main             = "Per-donor expression of consistent neutrophil DEGs (z-scored)",
    filename         = "Fig10_donor_expression_heatmap.pdf",
    width            = 5.5,
    height           = max(6, length(top_genes) * 0.18)
  )
  cat("\u2705 Saved: Fig10_donor_expression_heatmap.pdf\n")
} else {
  cat("\u26a0 Not enough consistent genes/donor columns for Fig 10 heatmap.\n")
}

# ============================================================
# FIGURE 11 — Per-donor mean-expression dotplot for top DEGs
# Each gene = a row; points = the 4 donors, coloured by condition.
# Shows the signal is carried by both diabetic donors, not one.
# ============================================================

if (length(top_genes) >= 2 && length(donor_cols) >= 2) {
  
  # Order genes by diabetic-minus-healthy mean (largest separation on top)
  gene_order <- donor_means %>%
    filter(gene %in% top_genes) %>%
    mutate(
      mean_D = rowMeans(across(all_of(intersect(DIABETIC_TAGS, donor_cols)))),
      mean_H = rowMeans(across(all_of(intersect(HEALTHY_TAGS,  donor_cols)))),
      delta  = mean_D - mean_H
    ) %>%
    arrange(delta) %>%
    pull(gene)
  
  # Limit to top 30 by |delta| for legibility
  dot_genes <- donor_means %>%
    filter(gene %in% top_genes) %>%
    mutate(
      mean_D = rowMeans(across(all_of(intersect(DIABETIC_TAGS, donor_cols)))),
      mean_H = rowMeans(across(all_of(intersect(HEALTHY_TAGS,  donor_cols)))),
      adelta = abs(mean_D - mean_H)
    ) %>%
    arrange(desc(adelta)) %>%
    slice_head(n = 30) %>%
    pull(gene)
  
  dot_df <- donor_means %>%
    filter(gene %in% dot_genes) %>%
    pivot_longer(all_of(donor_cols),
                 names_to = "donor", values_to = "mean_expr") %>%
    mutate(
      Condition = ifelse(donor %in% DIABETIC_TAGS, "Diabetic", "Healthy"),
      gene = factor(gene, levels = intersect(gene_order, dot_genes))
    )
  
  p11 <- ggplot(dot_df,
                aes(x = mean_expr, y = gene, colour = Condition)) +
    geom_line(aes(group = gene), colour = "grey80", linewidth = 0.5) +
    geom_point(aes(shape = donor), size = 2.6, alpha = 0.9) +
    scale_colour_manual(values = c(Healthy = "#4A90D9", Diabetic = NEUTRO_COLOR)) +
    labs(
      title    = "Per-Donor Mean Expression of Top Neutrophil DEGs",
      subtitle = "Each point = one donor (2 healthy, 2 diabetic) | line links the 4 donor means per gene",
      x = "Mean SCT expression (per donor)", y = NULL,
      colour = NULL, shape = "Donor"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, colour = "grey40"),
      axis.text.y   = element_text(size = 7, face = "italic"),
      legend.position = "right"
    )
  
  print(p11)
  ggsave("Fig11_donor_topDEG_dotplot.pdf", p11,
         width = 8, height = max(6, length(dot_genes) * 0.22), dpi = 300)
  cat("\u2705 Saved: Fig11_donor_topDEG_dotplot.pdf\n")
} else {
  cat("\u26a0 Not enough genes for Fig 11 dotplot.\n")
}

# ============================================================
# FIGURE 12 — DEGs recovered per named donor pairing
# x-axis = the actual 4 pairings (Diabetic donor vs Healthy donor)
# ============================================================
# Reads the per-pairing count CSV written by 07b.

if (file.exists("neutrophil_DEG_per_pairing_counts.csv")) {
  
  pc <- read.csv("neutrophil_DEG_per_pairing_counts.csv", stringsAsFactors = FALSE)
  
  # Friendly pairing labels: "D01 vs H04" style from the raw "..._vs_..." string
  short_tag <- function(x) {
    x <- gsub("SampleTag", "", x)
    x <- gsub("_hs", "", x)
    ifelse(paste0("SampleTag", x, "_hs") %in% DIABETIC_TAGS,
           paste0("D", x), paste0("H", x))
  }
  pc <- pc %>%
    tidyr::separate(pairing, into = c("diab_raw", "heal_raw"),
                    sep = "_vs_", remove = FALSE) %>%
    mutate(pair_label = paste0(short_tag(diab_raw), " vs ", short_tag(heal_raw)))
  
  # Long form for stacked Up/Down bars
  pc_long <- pc %>%
    select(pair_label, Up, Down) %>%
    tidyr::pivot_longer(c(Up, Down),
                        names_to = "Direction", values_to = "Count") %>%
    mutate(Direction = factor(Direction, levels = c("Up", "Down")),
           pair_label = factor(pair_label,
                               levels = pc$pair_label[order(-pc$Total)]))
  
  p12 <- ggplot(pc_long, aes(x = pair_label, y = Count, fill = Direction)) +
    geom_col(width = 0.65) +
    geom_text(data = pc,
              aes(x = factor(pair_label, levels = pc$pair_label[order(-pc$Total)]),
                  y = Total, label = Total),
              inherit.aes = FALSE, vjust = -0.4, size = 3.8,
              fontface = "bold", colour = "grey20") +
    scale_fill_manual(
      values = c("Up" = "#E05A5A", "Down" = "#4A90D9"),
      labels = c("Up" = "\u2191 Up in Diabetic", "Down" = "\u2193 Down in Diabetic")
    ) +
    expand_limits(y = max(pc$Total) * 1.12) +
    labs(
      title    = "DEGs Recovered per Donor Pairing",
      subtitle = "Each bar = one Diabetic-donor vs Healthy-donor comparison (D = diabetic, H = healthy)",
      x = NULL, y = "Number of significant DEGs (adj p < 0.05)", fill = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle   = element_text(size = 9.5, hjust = 0.5, colour = "grey40"),
      axis.text.x     = element_text(size = 11),
      panel.grid.major.y = element_line(colour = "grey90", linetype = "dashed")
    )
  
  print(p12)
  ggsave("Fig12_DEGs_per_pairing.pdf", p12,
         width = 7.5, height = 5.5, dpi = 300)
  cat("\u2705 Saved: Fig12_DEGs_per_pairing.pdf\n")
} else {
  cat("\u26a0 neutrophil_DEG_per_pairing_counts.csv not found \u2014 run 07b first.\n")
}

# ============================================================
# SUMMARY
# ============================================================
cat("  PER-DONOR VISUALISATION (08b) COMPLETE\n")
cat(sprintf("DEGs in consistency table : %d\n", nrow(consist)))
cat(sprintf("Recovered in all %d pairings (strict): %d\n",
            n_pairings, sum(consist$n_same_dir_sig == n_pairings)))
cat(sprintf("Direction-consistent (>= %d/%d sign agreement): %d\n",
            max(1, n_pairings - 1), n_pairings,
            sum(consist$n_same_sign >= max(1, n_pairings - 1))))
cat("Figures: Fig9, Fig10, Fig11\n")

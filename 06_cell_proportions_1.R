# ============================================================
# SECTION 06: CELL PROPORTION ANALYSIS
# ============================================================
# INPUT  : merged_annotated.rds
# OUTPUT : prop$data  (proportion table)
# SAVES  : cell_type_proportions.csv
# ============================================================

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(scMEGA)

# ── Load checkpoint ────────────────────────────────────────────────────────
merged <- readRDS("merged_annotated.rds")
cat("Loaded: merged_annotated.rds |", ncol(merged), "cells\n")

# Rebuild propcolors if not in environment
nb_cols    <- length(unique(na.omit(merged$Cell_annot)))
propcolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb_cols)
names(propcolors) <- sort(unique(na.omit(merged$Cell_annot)))

# ============================================================
# 6.1  PROPORTION BAR PLOT (scMEGA)
# ============================================================

prop <- CellPropPlot(merged,
                     group.by = "Cell_annot",
                     prop.in  = "Condition",
                     cols     = propcolors)
print(prop)
print(prop$data, n = 50)

# ============================================================
# 6.2  DONUT PLOT
# ============================================================

prop_donut <- ggplot(
  prop$data,
  aes(x = 3, y = proportion, fill = Cell_annot)
) +
  geom_col(width = 1.5, color = "white") +
  facet_grid(. ~ Condition) +
  coord_polar(theta = "y") +
  xlim(c(0.2, 3.8)) +
  scale_fill_manual(values = propcolors) +
  theme_void() +
  theme(
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(size = 15, face = "bold", hjust = 0.5)
  ) +
  geom_text(
    aes(label = ifelse(proportion > 0.05,
                       paste0(round(proportion * 100, 1), "%"), "")),
    position = position_stack(vjust = 0.6),
    size     = 3.5,
    color    = "black",
    fontface = "bold"
  ) +
  ggtitle("Cell Type Proportions — Healthy vs Diabetic")

print(prop_donut)
ggsave("Fig1_UMAP_prop.pdf", prop_donut,
       width = 10, height = 5.5, dpi = 300)

# ============================================================
# 6.3  EXPORT
# ============================================================

write.csv(prop$data, "cell_type_proportions.csv", row.names = FALSE)
cat("✅ Saved: cell_type_proportions.csv\n")
cat("✅ Section 06 complete\n")

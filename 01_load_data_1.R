# ============================================================
# SECTION 01: LOAD DATA
# scRNA-seq BD Rhapsody — Cart3 + Cart4
# Human Peripheral WBCs: Healthy vs Diabetic
# ============================================================
# INPUT  : Raw MEX directories + Sample Tag CSVs
# OUTPUT : merged  (Seurat object, tag-filtered, condition-labelled)
# SAVES  : (none — object passed to Section 02)
# ============================================================

library(Seurat)
library(dplyr)

set.seed(42)
options(future.globals.maxSize = 8000 * 1024^2)

# ── Paths ─────────────────────────────────────────────────────────────────
PATH_C3_COUNTS <- "C:/IISc/ScRNA_Analysis_Human/Cart3/CART3HUMAN_RSEC_MolsPerCell_MEX/"
PATH_C3_TAGS   <- "C:/IISc/ScRNA_Analysis_Human/Cart3/CART3HUMAN_Sample_Tag_Calls.csv"
PATH_C4_COUNTS <- "C:/IISc/ScRNA_Analysis_Human/Cart4/CART4HUMAN_RSEC_MolsPerCell_MEX/"
PATH_C4_TAGS   <- "C:/IISc/ScRNA_Analysis_Human/Cart4/CART4HUMAN_Sample_Tag_Calls.csv"

VALID_TAGS <- c(
  "SampleTag01_hs",   # Diabetic — Cart3
  "SampleTag04_hs",   # Healthy  — Cart3
  "SampleTag05_hs",   # Diabetic — Cart3
  "SampleTag12_hs"    # Healthy  — Cart4
)

# ============================================================
# 1.1  LOAD CART3
# ============================================================

counts_c3 <- Read10X(PATH_C3_COUNTS)

cart3 <- CreateSeuratObject(
  counts  = counts_c3[["Gene Expression"]],
  project = "Cart3"
)
cart3[["Protein"]] <- CreateAssayObject(
  counts = counts_c3$`Antibody Capture`
)

meta_c3 <- read.table(PATH_C3_TAGS,
                      sep = ",", header = TRUE, row.names = 1)
cart3 <- AddMetaData(cart3, metadata = meta_c3)

# Tag cartridge identity BEFORE merge — critical for Harmony later
cart3$Cartridge  <- "Cart3"
cart3$orig.ident <- "Cart3"

cat("Cart3 raw cells:", ncol(cart3), "\n")
print(table(cart3$Sample_Tag))

# ============================================================
# 1.2  LOAD CART4
# ============================================================

counts_c4 <- Read10X(PATH_C4_COUNTS)

cart4 <- CreateSeuratObject(
  counts  = counts_c4[["Gene Expression"]],
  project = "Cart4"
)
cart4[["Protein"]] <- CreateAssayObject(
  counts = counts_c4$`Antibody Capture`
)

meta_c4 <- read.table(PATH_C4_TAGS,
                      sep = ",", header = TRUE, row.names = 1)
cart4 <- AddMetaData(cart4, metadata = meta_c4)

cart4$Cartridge  <- "Cart4"
cart4$orig.ident <- "Cart4"

cat("Cart4 raw cells:", ncol(cart4), "\n")
print(table(cart4$Sample_Tag))

# ============================================================
# 1.3  MERGE
# add.cell.ids prevents barcode collision between cartridges
# ============================================================

merged <- merge(
  cart3,
  y            = cart4,
  add.cell.ids = c("C3", "C4"),
  project      = "BD_Rhapsody_WBC"
)

# merge() resets orig.ident — restore from Cartridge column
merged$orig.ident <- factor(as.character(merged$Cartridge))

cat("Any duplicate barcodes:", any(duplicated(colnames(merged))), "\n")
cat("\n=== All Sample Tags (pre-filter) ===\n")
print(table(merged$Sample_Tag))

# ============================================================
# 1.4  FILTER TO VALID SAMPLE TAGS ONLY
# Removes: Undetermined, Multiplet, unused tags
# ============================================================

merged <- subset(merged,
                 subset = Sample_Tag %in% VALID_TAGS)

cat("\n=== Sample Tags after filter ===\n")
print(table(merged$Sample_Tag))
cat("Total cells after tag filter:", ncol(merged), "\n")

# ============================================================
# 1.5  ADD CONDITION METADATA
# Using case_when() — robust to factor levels, no replace() risk
# ============================================================

merged$Condition <- dplyr::case_when(
  merged$Sample_Tag %in% c( "SampleTag01_hs", "SampleTag05_hs") ~ "Diabetic",
  merged$Sample_Tag %in% c( "SampleTag04_hs", "SampleTag12_hs") ~ "Healthy",
  TRUE ~ NA_character_
)
unique(merged$Condition)
merged$Condition <- factor(merged$Condition,
                           levels = c("Healthy", "Diabetic"))

cat("\n=== Condition Assignment ===\n")
cat("NA in Condition:", sum(is.na(merged$Condition)), "\n\n")
print(table(merged$Condition))
cat("\n=== Condition × Sample_Tag ===\n")
print(table(merged$Sample_Tag, merged$Condition))
cat("\n=== Condition × Cartridge ===\n")
print(table(merged$Cartridge, merged$Condition))

saveRDS(merged, "Merged.rds")
cat("\n✅ Section 01 complete — 'merged' object ready\n")

# ============================================================
# RNA-seq differential expression analysis with edgeR
# Time course: D0, D3, D6
# ============================================================

library(edgeR)
library(tidyverse)
library(pheatmap)
library(ggrepel)

# -------------------------
# 1. Define input files
# -------------------------

DATA_DIR <- "/home/course/rnaseq/data/"

COUNTS_FILE <- file.path(DATA_DIR, "counts/joint_withpairs_featureCounts.txt")
GENE_ANNOT_FILE <- file.path(DATA_DIR, "geneid2genename.txt")

OUTDIR <- "~/rnaseq/"

dir.create(OUTDIR, showWarnings = FALSE)

# -------------------------
# 2. Load count matrix
# -------------------------

# featureCounts files contain comment lines starting with "#".
# comment.char = "#" tells R to skip those lines.
fc <- read.delim(
  COUNTS_FILE,
  comment.char = "#",
  check.names = FALSE
)

# Select only count columns.
# These are the columns ending with ".bam".
count_cols <- grep(".bam$", colnames(fc), value = TRUE)

counts <- fc[, count_cols]

# Use gene IDs as row names.
rownames(counts) <- fc$Geneid

# Clean sample names to make them easier to read.
colnames(counts) <- colnames(counts) |>
  str_replace(".mm39.bam", "") |>
  str_replace("RNAseq_NCC_", "")

# -------------------------
# 3. Create sample metadata
# -------------------------

sample_info <- data.frame(
  sample = colnames(counts),
  timepoint = case_when(
    str_detect(colnames(counts), "D0") ~ "D0",
    str_detect(colnames(counts), "D3") ~ "D3",
    str_detect(colnames(counts), "D6") ~ "D6"
  )
)

# Set D0 as the reference level.
sample_info$timepoint <- factor(
  sample_info$timepoint,
  levels = c("D0", "D3", "D6")
)

rownames(sample_info) <- sample_info$sample

# Check that metadata matches count matrix.
print(sample_info)


# -------------------------
# 5. Create edgeR object
# -------------------------

group <- sample_info$timepoint

y <- DGEList(
  counts = counts,
  group = group
)

# -------------------------
# 6. Filter low-expression genes
# -------------------------

# Genes with very low counts are usually not informative.
# Removing them improves statistical power.
keep <- filterByExpr(y, group = group)

y <- y[keep, , keep.lib.sizes = FALSE]

# -------------------------
# 7. Normalize libraries
# -------------------------

# TMM normalization corrects for library size/composition differences.
y <- calcNormFactors(y)

# -------------------------
# 8. PCA for sample exploration
# -------------------------

logCPM <- cpm(y, log = TRUE, prior.count = 2)

pca <- prcomp(t(logCPM), scale. = FALSE)

percent_var <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  1
)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  timepoint = group,
  sample = colnames(logCPM)
)

p_pca <- ggplot(
  pca_df,
  aes(PC1, PC2, color = timepoint)
) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = sample), size = 3) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)"),
    color = "Time point"
  )

p_pca

ggsave(
  file.path(OUTDIR, "PCA_edgeR.pdf"),
  p_pca,
  width = 6,
  height = 5
)

# -------------------------
# 9. Differential expression model
# -------------------------

# The design matrix describes the experiment:
# expression ~ timepoint
design <- model.matrix(~ group)

# Estimate biological variability between replicates.
y <- estimateDisp(y, design)

# Fit the model.
fit <- glmQLFit(y, design)

# -------------------------
# 10. Global test: genes changing across D0, D3, D6
# -------------------------

# coef = 2:3 tests whether D3 or D6 differs from D0.
# This identifies genes that change across the time course.
qlf_time <- glmQLFTest(
  fit,
  coef = 2:3
)

edgeR_all <- topTags(qlf_time, n = Inf)$table |>
  rownames_to_column("gene_id")

# Load gene annotation
gene_annot <- read.delim(
  GENE_ANNOT_FILE,
  header = FALSE,
  col.names = c("gene_id", "gene_name", "gene_biotype")
)

gene_annot <- gene_annot |>
  mutate(
    gene_category = ifelse(
      gene_biotype == "protein_coding",
      "Coding",
      "Non-coding"
    )
  )
  
#rename gene id to gene name
edgeR_all <- topTags(qlf_time, n = Inf)$table |>
  rownames_to_column("gene_id") |>
  left_join(gene_annot, by = "gene_id") |>
  mutate(
    gene_label = ifelse(is.na(gene_name), gene_id, gene_name),
    comparison = "All_timepoints",
    result = ifelse(FDR < 0.05, "Dynamic", "Not significant")
  )

write.csv(
  edgeR_all,
  file.path(OUTDIR, "edgeR_all_timepoints_annotated.csv"),
  row.names = FALSE
)

# -------------------------
# 11. Pairwise comparisons
# -------------------------

contrast_matrix <- makeContrasts(
  D3_vs_D0 = groupD3,
  D6_vs_D0 = groupD6,
  D6_vs_D3 = groupD6 - groupD3,
  levels = design
)

pairwise_results <- list()

for (contrast_name in colnames(contrast_matrix)) {
  
  qlf <- glmQLFTest(
    fit,
    contrast = contrast_matrix[, contrast_name]
  )
  
  result <- topTags(qlf, n = Inf)$table |>
    rownames_to_column("gene_id") |>
    left_join(gene_annot, by = "gene_id") |>
    mutate(
      gene_label = ifelse(is.na(gene_name), gene_id, gene_name),
      comparison = contrast_name,
      regulation = case_when(
        FDR < 0.05 & logFC > 1 ~ "Up",
        FDR < 0.05 & logFC < -1 ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  pairwise_results[[contrast_name]] <- result
  
  write.csv(
    result,
    file.path(OUTDIR, paste0("edgeR_", contrast_name, "_annotated.csv")),
    row.names = FALSE
  )
}

pairwise_all <- bind_rows(pairwise_results)

# -------------------------
# 12. Bar plot: number of up/down genes
# -------------------------

deg_counts <- pairwise_all |>
  filter(regulation %in% c("Up", "Down")) |>
  count(comparison, regulation, gene_category)

p_bar <- ggplot(
  deg_counts,
  aes(x = comparison, y = n, fill = regulation)
) +
  geom_col(position = "dodge") +
  facet_wrap(~ gene_category) +
  theme_classic(base_size = 14) +
  labs(
    x = "Comparison",
    y = "Number of differentially expressed genes",
    fill = "Regulation"
  )

p_bar

ggsave(
  file.path(OUTDIR, "DEG_counts_barplot.pdf"),
  p_bar,
  width = 8,
  height = 5
)

# -------------------------
# 13. Volcano plots
# -------------------------

for (contrast_name in names(pairwise_results)) {
  
  result <- pairwise_results[[contrast_name]]
  
  p_volcano <- ggplot(
    result,
    aes(x = logFC, y = -log10(FDR), color = regulation)
  ) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_classic(base_size = 14) +
    labs(
      title = contrast_name,
      x = "log2 fold change",
      y = "-log10(FDR)",
      color = "Regulation"
    )
  
  print(p_volcano)
  
  ggsave(
    file.path(OUTDIR, paste0("volcano_", contrast_name, ".pdf")),
    p_volcano,
    width = 6,
    height = 5
  )
}

# -------------------------
# 14. Heatmap of top dynamic genes
# -------------------------

top_genes <- edgeR_all |>
  filter(!is.na(FDR)) |>
  arrange(FDR) |>
  slice_head(n = 50) |>
  pull(gene_id)

mat <- logCPM[top_genes, ]

# Center each gene around its average expression.
# This makes patterns across samples easier to see.
mat <- mat - rowMeans(mat)

gene_labels <- edgeR_all |>
  filter(gene_id %in% top_genes) |>
  arrange(match(gene_id, top_genes)) |>
  pull(gene_label)

rownames(mat) <- make.unique(gene_labels)

plot_heat <- pheatmap(
  mat,
  annotation_col = sample_info[, "timepoint", drop = FALSE],
  show_rownames = TRUE,
  fontsize_row = 6
)

plot_heat

pdf(
  file.path(OUTDIR, "top50_dynamic_genes_heatmap.pdf"),
  width = 7,
  height = 9
)

print(plot_heat)

dev.off()

# -------------------------
# 15. GO enrichment of top dynamic genes
# -------------------------

# This analysis asks:
# What biological processes are enriched among the top dynamic genes?

library(clusterProfiler)
library(org.Mm.eg.db)

top_genes_clean <- sub("\\..*", "", top_genes)

# Convert Ensembl gene IDs to Entrez IDs.
# clusterProfiler uses Entrez IDs for GO enrichment.
gene_conversion <- bitr(
  top_genes_clean,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Mm.eg.db
)

# Run GO enrichment for Biological Process.
go_results <- enrichGO(
  gene = gene_conversion$ENTREZID,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save full GO table.
write.csv(
  as.data.frame(go_results),
  file.path(OUTDIR, "GO_top50_dynamic_genes.csv"),
  row.names = FALSE
)

# Dotplot of enriched GO terms.
p_go <- dotplot(
  go_results,
  showCategory = 15
) +
  ggtitle("GO enrichment: top dynamic genes")

p_go

ggsave(
  file.path(OUTDIR, "GO_top50_dynamic_genes_dotplot.pdf"),
  p_go,
  width = 8,
  height = 6
)
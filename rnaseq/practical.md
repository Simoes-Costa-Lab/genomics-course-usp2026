---
title: Practical Session
parent: Transcriptomics (RNA-seq)
nav_order: 2
---

# Practical Session

In this practical session, we will analyze a real RNA-seq dataset using the **edegR** package in R.

Our goal is to follow a typical differential RNA-seq workflow, starting from a count matrix and ending with the identification of genes that are differentially expressed between conditions.

# Learning Objectives

By the end of this session, participants will be able to:

- load and inspect a RNA-seq dataset
- perform basic quality control
- generate low-dimensional representations of samples
- normalize expression data
- identify differentially expressed genes

# About this Data

The dataset used in this exercise is a published dataset from this publication: [BRD4 binds to active cranial neural crest enhancers to regulate RUNX2 activity during osteoblast differentiation](https://pubmed.ncbi.nlm.nih.gov/38063851/) by Musa et al. in Development (2024). 

In this study, the authors characterized the cellular and molecular function of the gene BRD4 for craniofacial development. Of all the data they generated, we are particularly interested in the osteogenic differentiation of a cell line called O9-1. This is a cranial neural crest cell that can be differentiated into osteoblast. Since they conducted differentiation in wild type cell lines, we have 3 replicates of the most important time points.

	
|id|sample|timepoint|
|------|------|----------|
|D0_WT_1|D0_WT_1|D0|
|D0_WT_2|D0_WT_2|D0|
|D0_WT_3|D0_WT_3|D0|
|D3_OST_WT_1|D3_OST_WT_1|D3|
|D3_OST_WT_2|D3_OST_WT_2|D3|
|D3_OST_WT_3|D3_OST_WT_3|D3|
|D6_OST_WT_1|D6_OST_WT_1|D6|
|D6_OST_WT_2|D6_OST_WT_2|D6|
|D6_OST_WT_3|D6_OST_WT_3|D6|

# Workflow Overview

<div align="center">
<img src="/genomics-course-usp2026/assets/images/rnaseq_bioinfo.png" width="500">
</div>

<p align="center">
<em>Figure 1. Overview of the RNA-seq workflow. Source: https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/ </em>
</p>

# Step 1 — Loading Packages

```r
library(edgeR)
library(tidyverse)
library(pheatmap)
```

# Step 2 — Set input and output paths

```r
# Set path. It's useful to do it at the beginning, so you don't have to rewrite file paths again and again
DATA_DIR <- "/home/course/rnaseq/data/"

COUNTS_FILE <- file.path(DATA_DIR, "counts/joint_withpairs_featureCounts.txt")
GENE_ANNOT_FILE <- file.path(DATA_DIR, "geneid2genename.txt")

OUTDIR <- "~/rnaseq/"

dir.create(OUTDIR, showWarnings = FALSE)
```

# Step 3 — Loading Data

```r
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
```

# Step 4 — Make/Read metadata

Up until now we have loaded a gene expression matrix, where each column is a sample and each row is a gene. To be able to do a differential expression analysis, we need to create a metadata table, so that edegR can group samples of the same condition. Here, our sample name already tells us everything we need, so simple parsing using R is enough. However, I would recommend you to build a metadata csv file to be sure your sample id match the condition, and then just read it into your R script. 

```r
# Make metadata
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
```

Now we can check our metadata and save it for future reference:

```r
# Check
print(sample_info)

write.csv(
  sample_info,
  file.path(OUTDIR, "metadata.csv"),
  row.names = FALSE
)
```

# Step 3 — Differential Expression Analysis with edgeR

The goal of this analysis is to identify genes whose expression changes during osteoblast differentiation. Our experiment contains three developmental stages with three biological replicates each:

* **D0**: Neural crest cells
* **D3**: Early osteoblast differentiation
* **D6**: Late osteoblast differentiation

## Define Experimental Groups

First, we specify the experimental condition associated with each sample. This information will be used throughout the analysis to compare gene expression between time points.

```r
group <- sample_info$timepoint
```

## Create an edgeR Object

```r
y <- DGEList(
  counts = counts,
  group = group
)
```

The count matrix is stored inside a `DGEList` object.This object contains:

* the raw count matrix
* sample information
* normalization factors
* statistical parameters estimated during the analysis



## Filter Lowly Expressed Genes

Genes with extremely low expression provide little statistical power and can increase noise.

```r
keep <- filterByExpr(y, group = group)

y <- y[keep, , keep.lib.sizes = FALSE]
```

## Normalize Sequencing Depth

Different samples often have different library sizes. edgeR uses the **TMM (Trimmed Mean of M-values)** method to estimate normalization factors. This allows gene expression levels to be compared more fairly across samples.

```r
y <- calcNormFactors(y)
```

# Principal Component Analysis (PCA)

Now that we have normalized our data, this is a great point to stop and inspect our samples. Before doing complex calculations we can have a general sense of similarity/differences between samples using Principal Components Analysis: a dimensionality-reduction method employed to stratify data by their variance. In other words, points that are more closely together on a PCA are more similar to each other, while points that are more dissimilar to each other will be further apart. 

```r
# edgeR logCPM for visualization
logCPM <- cpm(y, log = TRUE, prior.count = 2)

# PCA
pca <- prcomp(t(logCPM), scale. = FALSE)
```

To visualize the PCA we have to create a new dataframe:

```r
# Percentage of variance explained
percent_var <- round(
  100 * (pca$sdev^2 / sum(pca$sdev^2)),
  1
)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  timepoint = group,
  sample = colnames(logCPM)
)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  timepoint = group,
  sample = colnames(logCPM)
)
```

We can visualize this using a scatter plot:

```r
p_pca <- ggplot(
  pca_df,
  aes(PC1, PC2, color = timepoint)
) +
  geom_point(size = 4) +
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
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot01.png" width="500">
</div>

<p align="center">
<em>Figure 2. Dimension plot of the first 2 PCs </em>
</p>

</details>

## Differential expression model

```r
# The design matrix describes the experiment:
# expression ~ timepoint
design <- model.matrix(~ group)

# Estimate biological variability between replicates.
y <- estimateDisp(y, design)

# Fit the model.
fit <- glmQLFit(y, design)
```

Biological replicates are never identical. edgeR estimates gene-specific variability using dispersion parameters. Dispersion reflects the amount of variation observed between replicates. Genes with higher dispersion show greater biological variability.

Next, edgeR fits a generalized linear model (GLM) to every gene. This model estimates how gene expression changes across the different time points.

## Test for Differential Expression

To identify genes that change during differentiation, we perform a quasi-likelihood F-test.

```r
# coef = 2:3 tests whether D3 or D6 differs from D0.
# This identifies genes that change across the time course.
qlf_time <- glmQLFTest(
  fit,
  coef = 2:3
)
```

This test asks:

> Does this gene show a significant change in expression at any point during differentiation?

Unlike a pairwise comparison, this test evaluates the entire trajectory:

```text
D0 → D3 → D6
```

Genes significant in this analysis are considered dynamically regulated during osteoblast differentiation.

## Extract Results

Finally, we retrieve all genes ranked by statistical significance.

```r
edgeR_all <- topTags(
  qlf_time,
  n = Inf
)$table |>
  rownames_to_column("gene_id")
```

The resulting table contains:

| Column | Description                        |
| ------ | ---------------------------------- |
| logFC  | Estimated fold change              |
| logCPM | Average expression level           |
| F      | Test statistic                     |
| PValue | Raw p-value                        |
| FDR    | Multiple-testing corrected p-value |

Genes with low FDR values are strong candidates for involvement in the differentiation process.

At this stage, investigate the table. Do you know any of those genes? Probably not, because we are using gene id up until this point. To make the upcoming analysis more interpretable for a human, let's rename the gene id to gene name using an auxiliary file:

```r
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
```

## Pairwise Comparisons

After identifying genes that change across the entire time course, we can perform specific pairwise comparisons.

```r
contrast_matrix <- makeContrasts(
  D3_vs_D0 = groupD3,
  D6_vs_D0 = groupD6,
  D6_vs_D3 = groupD6 - groupD3,
  levels = design
)
```

These contrasts allow us to ask more specific biological questions:

* Which genes change early during differentiation? (**D3 vs D0**)
* Which genes change late during differentiation? (**D6 vs D0**)
* Which genes continue changing between D3 and D6? (**D6 vs D3**)

To save the differentially expressed genes of each pairwise comparison, we can loop our contrast matrix and save the csv for future analysis.

```r
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
```

## Explore results through visualization

Let's extract some numbers: how many genes are up or down regulated in each pair-wise comparison? What type of genes are they?

```r
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
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot02.png" width="700">
</div>

<p align="center">
<em>Figure 3. Number of differentially expressed gene by gene type. </em>
</p>

</details>

One common way to visualize differentially expressed genes is using a type of scatterplot called volcano plot: where the x axis is the log fold change (negative and positive) and the y axis is the FDR value: the highest, the more you can trust the log fold change value. 

```r
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
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot03.png" width="500">
</div>

<p align="center">
<em>Figure 4. Volcano plot of Day 3 vs Day 0. </em>
</p>

</details>

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot04.png" width="500">
</div>

<p align="center">
<em>Figure 5. Volcano plot of Day 6 vs Day 0.  </em>
</p>

</details>

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot05.png" width="500">
</div>

<p align="center">
<em>Figure 6. Volcano plot of Day 6 vs Day3. </em>
</p>

</details>

Cool! We get to see all genes at once using the volcano plot. This is a good sanity check for your analysis, but is still very hard to assess a biological interpretation. We saw that we have more than  2000 genes either up or down regulated in each pairwise comparison. One way to extract biological meaning is looking at a slice of our data. Here we are choosing the top 50 most significantly differentially expressed gene in any pairwise comparison. 

```r
top_genes <- edgeR_all |>
  filter(!is.na(FDR)) |>
  arrange(FDR) |>
  slice_head(n = 50) |>
  pull(gene_id)
```

Now, let's build our gene matriz and calculate the z-score. With this we can plot a heatmap and see how genes are behaving across the time points.

```r
# Center each gene around its average expression.
# This makes patterns across samples easier to see.
mat <- logCPM[top_genes, ]

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
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot06.png" width="500">
</div>

<p align="center">
<em>Figure 7. Heatmap of top 50 most significantly differentially expressed gene in any pairwise comparison. </em>
</p>

</details>

This visualization is very useful, specially if you know your biology. We can see a lot of collagen genes showing up and the Ogn gene peaking at day 6. However, if you know nothing about this cell line, you probably saw the gene list and though: how is this any different from having the gene id? I still can't interpret it!

Another very useful analysis is a Gene Ontology enrichment of genes. Since we want to know more about our 50 genes, we will focus on them.

```r
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
```

We can also plot this result in a dotplot of GO terms (very nice for publications)

```r
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
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/rna/Rplot07.png" width="500">
</div>

<p align="center">
<em>Figure 7. GO terms of our top 50 differentially expressed genes. </em>
</p>

</details>

# The End

I hope this give you a nice base to think about your own analysis: we worked with a real dataset that has an experiment design of 3 time points.  As you can see, after doing a differential expression analysis we can ask a lot of questions, so I would encourage you to explore them on your own:
- What genes go down regulated from day 0 to day 6 in a sustained fashion? What about up regulated? 
- Is there a subset of genes that loose expression on day 3 but regain it on day 6? What about the other way around?
- What if I selected the top 50 genes with the greatest logFC value? Would my GO change? Why? And can I trust that?
- How can I use this analysis to generate further hypothesis for my research?

Please go to the  [References](./references)  tab for more resources!

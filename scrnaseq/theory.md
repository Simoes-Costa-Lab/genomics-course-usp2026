---
title: Theory
parent: Single‑cell RNA-seq
nav_order: 1
---

# Theory

# What is scRNA-seq?

Single-cell RNA-seq (scRNA-seq) is a technique used to measure gene expression individually in thousands of cells simultaneously.

Unlike conventional RNA-seq (bulk RNA-seq), where the observed signal represents the average expression across a cell population, scRNA-seq enables the investigation of cellular heterogeneity.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrnaseq_overview.png" width="700">
</div>

<p align="center">
<em>Figure 1. Overview of the single-cell RNA-seq workflow. Source: https://doi.org/10.1371/journal.ppat.1011898 </em>
</p>

---

# Common Biological Questions

scRNA-seq can be used to investigate:

- tissue cellular composition
- embryonic development
- tumor heterogeneity
- cellular responses to stimuli
- differentiation trajectories
- transient cellular states

---

# Bulk RNA-seq vs scRNA-seq

In bulk RNA-seq, the measured signal represents the average expression across thousands or millions of cells.

In contrast, scRNA-seq measures gene expression at the level of individual cells.

## Comparison

| Bulk RNA-seq | scRNA-seq |
|---|---|
| population average | single-cell resolution |
| lower technical noise | higher variability |
| lower cost | higher computational cost |
| average expression | cellular heterogeneity |

<div align="center">
<img src="/genomics-course-usp2026/assets/images/bulk_vs_singlecell.png" width="700">
</div>

<p align="center">
<em>Figure 2. Comparison between bulk RNA-seq and scRNA-seq. Source: https://www.completegenomics.com/methods/single-cell-rna-sequencing/ </em>
</p>

---

# Cell Isolation and Barcoding

In platforms such as 10x Genomics, individual cells are encapsulated into droplets containing beads with molecular barcodes.

Each RNA molecule receives:

- a cell barcode
- a UMI (Unique Molecular Identifier)

## Cell Barcodes

Cell barcodes allow identification of the cell from which each read originated.

## UMIs

UMIs enable the identification of unique RNA molecules and help reduce PCR amplification bias.

---

# Overview of the scRNA-seq Workflow

## General Workflow

<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrnaseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figure 3. Overview of a typical single-cell RNA-seq workflow.</em>
</p>

---

# Cell-by-Gene Matrix

After alignment and quantification, the data are organized into a matrix:

| Gene | Cell_1 | Cell_2 | ... | Cell_n |
|---|---|---|---|---|
| SOX10 | 4 | 0 | ... | 2 |
| PAX7 | 12 | 3 | ... | 0 |
| TFAP2A | 0 | 8 | ... | 23 |

In this matrix:

- rows represent genes
- columns represent cells
- values represent RNA abundance

---

# Sparsity and Dropouts

scRNA-seq matrices are typically highly sparse.

Many values appear as zero because of:

- low RNA capture efficiency
- limited sequencing depth
- very low expression levels
- technical limitations

These zero values are commonly referred to as **dropouts**.

---

# Quality Control

Quality control in scRNA-seq is performed at the level of individual cells.

The goal is to remove low-quality cells before downstream analysis.

## Common Metrics

| Metric | Interpretation |
|---|---|
| number of detected genes | cellular complexity |
| total number of UMIs | sequencing depth |
| percentage of mitochondrial transcripts | cellular integrity |

---

# Low-Quality Cells

Problematic cells frequently exhibit:

- few detected genes
- low total counts
- high percentages of mitochondrial transcripts

These cells often correspond to:

- dead cells
- debris
- empty droplets

---

# Doublets

In some cases, two cells may be encapsulated within the same droplet.

These events are known as **doublets**.

Doublets can generate artificial expression profiles that combine signals from two different cell types.

Common tools include:

- DoubletFinder
- Scrublet

---

# Normalization

Individual cells are sequenced at different depths.

Normalization aims to correct these technical differences before comparing cells.

## Goals of Normalization

- correct sequencing depth differences
- stabilize variance
- enable comparisons across cells

## Common Methods

| Method | Tool |
|---|---|
| LogNormalize | Seurat |
| SCTransform | Seurat |
| CPM normalization | classical approaches |

---

# Selection of Highly Variable Genes

Not all genes are equally informative.

Highly variable genes help identify:

- cellular differences
- biological states
- distinct populations

These genes are often used in downstream analyses.

---

# Dimensionality Reduction

scRNA-seq experiments measure thousands of genes per cell.

Dimensionality reduction techniques help summarize these high-dimensional datasets.

## Common Methods

| Method | Purpose |
|---|---|
| PCA | linear dimensionality reduction |
| UMAP | visualization |
| t-SNE | local structure visualization |

<div align="center">
<img src="/genomics-course-usp2026/assets/images/umap_example.png" width="700">
</div>

<p align="center">
<em>Figure 4. Example of a UMAP projection from an scRNA-seq experiment. Source: https://www.nature.com/articles/s41586-023-05869-0 </em>
</p>

---

# Cell Clustering

After dimensionality reduction, similar cells can be grouped together.

These clusters often correspond to:

- cell types
- cell states
- distinct biological populations

---

# Cell Annotation

After identifying clusters, each population must be interpreted biologically.

This is typically achieved using:

- known marker genes
- reference datasets
- published literature

---

# Differential Expression Analysis

scRNA-seq also enables comparisons of gene expression between:

- clusters
- experimental conditions
- cellular states

---

# Limitations of scRNA-seq

- high sparsity
- technical dropouts
- increased experimental noise
- substantial computational requirements
- loss of spatial information
- cell dissociation may alter gene expression
---
title: Theory
parent: ATAC-seq
nav_order: 1
---

# Theory

# What is ATAC-seq?

ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) is a technique used to identify accessible chromatin regions within a cell population.

The method uses a hyperactive transposase (Tn5) that preferentially fragments open chromatin regions while simultaneously inserting sequencing adapters.

Accessible regions identified by ATAC-seq frequently correspond to:

* enhancers
* promoters
* active regulatory regions
* transcription factor binding sites

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_overview.png" width="700">
</div>

<p align="center">
<em>Figure 1. Overview of the ATAC-seq technique. BioRender template.</em>
</p>

---

# Common Biological Questions

ATAC-seq can be used to investigate:

* active regulatory regions
* epigenetic responses to stimuli
* regulatory differences between tissues
* potentially active transcription factors

---

# Open and Closed Chromatin

Nuclear DNA is organized into structures called nucleosomes.

Highly compacted regions tend to be less accessible to the transcriptional machinery.

In contrast, open chromatin regions allow:

* transcription factor binding
* RNA polymerase recruitment
* gene activation

ATAC-seq exploits these differences in chromatin accessibility.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/chromatin_accessibility.png" width="700">
</div>

<p align="center">
<em>Figure 2. Open and closed chromatin regions. Source: https://www.nature.com/articles/s41576-018-0089-8</em>
</p>

---

# Experimental Design

As with RNA-seq, experimental design is critical for ATAC-seq experiments.

### Important Concepts

* biological replicates
* nuclear quality
* number of cells
* batch effects
* sequencing depth

## Important Considerations

ATAC-seq is particularly sensitive to:

* cellular degradation
* excessive lysis
* mitochondrial contamination
* poor nuclear quality

---

# Characteristics of ATAC-seq Data

In ATAC-seq:

* short fragments typically correspond to nucleosome-free regions
* larger fragments may reflect nucleosome positioning
* mitochondrial reads are often abundant

## Fragment Size Distribution

One of the most characteristic quality metrics of ATAC-seq.

Small fragments generally correspond to open chromatin regions free of nucleosomes.

Larger fragments may reflect:

* mono-nucleosomes
* di-nucleosomes
* higher-order chromatin organization

<div align="center">
<img src="/genomics-course-usp2026/assets/images/fragment_distribution.png" width="700">
</div>

<p align="center">
<em>Figure 3. Fragment size distribution in ATAC-seq. Source: https://www.activemotif.com/blog-library-qc</em>
</p>

---

# Overview of the ATAC-seq Workflow

A typical ATAC-seq experiment involves several computational steps.

## General Workflow

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figure 4. Overview of a typical ATAC-seq workflow.</em>
</p>

---

# Reads and FASTQ Files

As in RNA-seq, sequencing generates FASTQ files containing DNA reads.

Each read represents an accessible DNA fragment identified by the Tn5 transposase.

---

# Quality Control

Quality control is essential for evaluating:

* enrichment in open chromatin regions
* library quality
* signal-to-noise ratio

Common tools include:

* FastQC
* MultiQC
* deepTools

---

# Important ATAC-seq Metrics

## Read Quality

Evaluates:

* per-base quality
* GC content
* adapter contamination
* duplication levels

---

## Alignment Rate

Measures how many reads align correctly to the reference genome.

Low alignment rates may indicate:

* contamination
* poor sequencing quality
* library preparation issues

---

## Mitochondrial Reads

ATAC-seq frequently generates a large proportion of reads derived from mitochondrial DNA.

High proportions of mitochondrial reads may indicate:

* excessive cell lysis
* poor nuclear quality

---

# Alignment

After quality control, reads are aligned to a reference genome.

Common tools include:

| Tool    | Characteristics           |
| ------- | ------------------------- |
| Bowtie2 | Widely used for ATAC-seq  |
| BWA     | Efficient for short reads |

The goal is to determine the genomic position of accessible fragments.

## BAM Files

After alignment, reads are stored in BAM files containing:

* genomic coordinates
* orientation
* alignment quality
* paired-end information

---

# Duplicate Removal

During PCR amplification, fragments may be over-amplified.

Because ATAC-seq analyzes DNA rather than RNA, duplicated fragments are often technical artifacts rather than biological signals.

Duplicate reads can artificially inflate the signal.

Common tools include:

* samtools
* Picard

---

# Peak Calling

One of the central steps of ATAC-seq analysis.

The goal of peak calling is to identify genomic regions with significant enrichment of sequencing reads.

These regions represent potentially accessible regulatory elements.

## Common Tools

| Tool        | Characteristics           |
| ----------- | ------------------------- |
| MACS2/MACS3 | Most widely used standard |
| Genrich     | Optimized for ATAC-seq    |

## Peak Concept

A peak represents:

```text
high local fragment density
→ increased chromatin accessibility
```

<div align="center">
<img src="/genomics-course-usp2026/assets/images/peak_calling.png" width="700">
</div>

<p align="center">
<em>Figure 5. Example of peak identification in ATAC-seq. Source: https://www.nature.com/articles/s41467-025-67491-0</em>
</p>

---

# Accessibility Quantification

After peaks have been identified, chromatin accessibility can be quantified across samples.

The result is typically a matrix:

| Peak         | Sample_1 | Sample_2 |
| ------------ | -------- | -------- |
| chr1:1-100   | 120      | 340      |
| chr2:200-300 | 540      | 210      |

In this matrix:

* rows represent accessible regions
* columns represent samples
* values represent fragment abundance

---

# Normalization

As in RNA-seq, raw ATAC-seq counts are not directly comparable between samples.

Technical differences may arise due to:

* sequencing depth
* transposition efficiency
* library composition
* proportion of mitochondrial reads
* total number of detected peaks

Normalization aims to correct these differences and enable more reliable biological comparisons.

## Common Approaches

| Method              | Purpose                              |
| ------------------- | ------------------------------------ |
| CPM                 | Correct sequencing depth             |
| TMM                 | Robust normalization between samples |
| DESeq2 size factors | Statistical modeling of count data   |

After normalization, accessibility matrices can be used for downstream analyses such as:

* PCA
* clustering
* differential accessibility analysis
* heatmaps

---

# Differential Accessibility

ATAC-seq can be used to identify differentially accessible regions between conditions.

These analyses allow researchers to investigate:

* regulatory activation
* chromatin remodeling
* epigenetic changes

Common tools include:

* DESeq2
* edgeR
* DiffBind

---

# Limitations of ATAC-seq

* does not directly measure protein binding
* accessibility does not necessarily imply functional activity
* repetitive regions can complicate alignment
* highly sensitive to sample quality
* limited resolution in heterogeneous cell populations

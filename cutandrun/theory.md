---
title: Theory
parent: CUT&RUN
nav_order: 1
---

# Theory

# What is CUT&RUN?

CUT&RUN (Cleavage Under Targets and Release Using Nuclease) is a technique used to map protein–DNA interactions across the genome.

The approach uses specific antibodies to direct a nuclease to DNA regions associated with proteins of interest, such as:

* transcription factors
* modified histones
* regulatory proteins

After DNA cleavage at the target regions, the released fragments are sequenced and analyzed computationally.

CUT&RUN allows researchers to identify:

* transcription factor binding sites
* epigenetic marks
* active or repressed regulatory regions

<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutandrun_overview.png" width="500">
</div>

<p align="center">
<em>Figure 1. Overview of the CUT&RUN technique. Source: https://elifesciences.org/articles/21856</em>
</p>

---

# Common Biological Questions

CUT&RUN can be used to investigate:

* transcription factor binding sites
* distribution of histone modifications
* active regulatory elements
* mechanisms of gene regulation
* epigenetic remodeling

---

# Regulatory Proteins and Chromatin

Gene expression is controlled by proteins capable of interacting with specific DNA regions.

These proteins include:

* transcription factors
* cofactors
* chromatin remodeling proteins
* modified histones

CUT&RUN enables direct mapping of these interactions across the genome.

---

# Histone Modifications

Histones can undergo different chemical modifications associated with distinct regulatory states.

## Common Examples

| Mark     | Biological Association |
| -------- | ---------------------- |
| H3K27ac  | Active enhancers       |
| H3K4me3  | Active promoters       |
| H3K27me3 | Gene repression        |
| H3K4me1  | Poised enhancers       |

---

# Experimental Design

As with RNA-seq and ATAC-seq, experimental design is fundamental.

## Important Concepts

* biological replicates
* antibody quality
* negative controls
* sequencing depth
* batch effects

## Important Controls

CUT&RUN experiments frequently include:

* IgG control
* no-antibody control
* input DNA (less common)

---

# Overview of the CUT&RUN Workflow

A typical CUT&RUN experiment involves several computational steps.

### General Workflow

<div align="center">
<img src="/genomics-course-usp2026/assets/images/atacseq_workflow.png" width="700">
</div>

<p align="center">
<em>Figure 2. Overview of a typical CUT&RUN workflow.</em>
</p>

---

# Reads and FASTQ Files

Sequencing generates FASTQ files containing reads derived from genomic regions bound by the protein of interest.

In CUT&RUN:

* fragments are typically short
* background signal is generally low
* enrichment tends to be highly specific

---

# Quality Control

Quality control in CUT&RUN evaluates:

* read quality
* enrichment efficiency
* experimental background
* library complexity

Common tools include:

* FastQC
* MultiQC
* deepTools

---

# Important Metrics

## Read Quality

Evaluates:

* per-base quality
* GC content
* adapter contamination
* duplication levels

---

## Alignment Rate

Indicates how many reads align correctly to the reference genome.

Low alignment rates may indicate:

* contamination
* poor library quality
* excessively short fragments

---

## Fragment Size Distribution

CUT&RUN generally produces smaller and more specific fragments than ChIP-seq.

Abnormal fragment distributions may suggest:

* over-digestion
* poor experimental efficiency
* DNA degradation

---

# Alignment

After quality control, reads are aligned to the reference genome.

Common tools include:

| Tool    | Characteristics           |
| ------- | ------------------------- |
| Bowtie2 | Widely used               |
| BWA     | Efficient for short reads |

The goal is to identify where the fragments originated in the genome.

---

# BAM Files

After alignment, reads are stored in BAM files containing:

* genomic coordinates
* orientation
* alignment quality
* paired-end information

These files can be visualized using genome browsers such as:

* IGV
* UCSC Genome Browser

---

# Peak Calling

Peak calling is the step used to identify regions enriched for sequencing reads.

These regions represent potential binding sites for the protein being analyzed.

### Common Tools

| Tool        | Characteristics       |
| ----------- | --------------------- |
| MACS2/MACS3 | Widely used           |
| SEACR       | Optimized for CUT&RUN |

---

# Peak Concept

A peak represents:

```text
high local fragment density
→ potential protein–DNA interaction
```

<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutandrun_peaks.png" width="700">
</div>

<p align="center">
<em>Figure 3. Example of peaks identified in CUT&RUN. Source: https://genome.cshlp.org/content/30/1/35</em>
</p>

---

# Narrow Peaks and Broad Peaks

Different proteins generate distinct enrichment patterns.

### Narrow Peaks

Typically associated with:

* transcription factors
* localized binding events

## Broad Peaks

Typically associated with:

* histone modifications
* large chromatin domains

---

# Signal Quantification

After peaks have been identified, signal intensity can be quantified across samples.

The result is typically a matrix:

| Peak         | Sample_1 | Sample_2 |
| ------------ | -------- | -------- |
| chr1:1-100   | 120      | 340      |
| chr2:200-300 | 540      | 210      |

In this matrix:

* rows represent enriched regions
* columns represent samples
* values represent fragment abundance

---

# Normalization

In CUT&RUN, technical differences between libraries can influence the observed signal intensity across samples.

These differences include:

* sequencing depth
* digestion efficiency
* antibody enrichment efficiency
* library complexity
* experimental background

Normalization aims to reduce these technical effects before biological comparisons.

## Common Approaches

| Method              | Purpose                          |
| ------------------- | -------------------------------- |
| CPM                 | Correct sequencing depth         |
| RPGC                | Normalize by genomic coverage    |
| DESeq2 size factors | Robust comparison across samples |

After normalization, the data can be used for:

* enrichment comparisons
* heatmap generation
* metaplot visualization
* differential peak analysis

---

# Differential Enrichment

CUT&RUN can be used to compare enrichment patterns between biological conditions.

These analyses allow researchers to investigate:

* gain or loss of binding
* epigenetic remodeling
* regulatory changes

Common tools include:

* DiffBind
* DESeq2
* edgeR

---

# Motif Analysis

After identifying peaks, enriched DNA sequence motifs can be identified.

Motifs may indicate:

* active transcription factors
* regulatory programs
* potential gene regulatory networks

Common tools include:

* HOMER
* MEME
* chromVAR

---

# Integration with RNA-seq and ATAC-seq

CUT&RUN is frequently integrated with other genomic approaches.

## Examples

| Technique | Information             |
| --------- | ----------------------- |
| RNA-seq   | Gene expression         |
| ATAC-seq  | Chromatin accessibility |
| CUT&RUN   | Protein–DNA binding     |

Integrating these layers allows researchers to investigate regulatory mechanisms more comprehensively.

---

# Limitations of CUT&RUN

* dependence on antibody quality
* resolution depends on the protein being analyzed
* experimental background may still occur
* low-abundance proteins may generate weak signals
* functional interpretation is not always straightforward

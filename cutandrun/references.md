---
title: References
parent: CUT&RUN
nav_order: 3
---

# References

This page gathers recommended datasets, tutorials, software documentation, and additional resources for learning and exploring CUT&RUN analysis.

---

# Dataset Used in This Practical Session

---

# Scripts for this Module

---

# Tutorials

### CUT&RUN Complete Guide

Overview of the method, experimental design, advantages, limitations, and analysis strategies.

- https://www.activemotif.com/documents/CUT-RUN-Complete-Guide.pdf

### Choosing Genomics Tools: CUT&RUN and CUT&Tag

Excellent overview comparing CUT&RUN, CUT&Tag, and ChIP-seq.

- https://hutchdatascience.org/Choosing_Genomics_Tools/cutrun-and-cuttag.html

---

# Software Documentation

### MACS3

Peak calling software commonly used for CUT&RUN and ATAC-seq.

- https://macs3-project.github.io/MACS/

### SEACR

Peak caller specifically developed for CUT&RUN datasets.

- http://seacr.fredhutch.org

### HOMER

Motif discovery and peak annotation.

- http://homer.ucsd.edu/homer

### ChIPseeker

Peak annotation and visualization.

- https://bioconductor.org/packages/ChIPseeker

### deepTools

Signal visualization and QC.

- https://deeptools.readthedocs.io

### IGV

Genome browser for visualizing genomic tracks.

- https://igv.org

---

# Further Reading

## CUT&RUN Methodology

### An Efficient Targeted Nuclease Strategy for High-Resolution Mapping of DNA Binding Sites

The original CUT&RUN publication introducing the method.

Skene & Henikoff, 2017.

- https://elifesciences.org/articles/21856

### Improved CUT&RUN Chromatin Profiling Tools

Protocol improvements and optimization strategies.

Meers et al., 2019.

- https://elifesciences.org/articles/46314

---

## Peak Calling

### Peak Calling by Sparse Enrichment Analysis for CUT&RUN Chromatin Profiling

Introduces SEACR, a peak caller designed specifically for low-background CUT&RUN datasets.

- https://doi.org/10.1186/s13072-019-0287-4

### Benchmarking Peak Calling Methods for CUT&RUN

Comparison of SEACR, MACS2, GoPeaks and other peak callers.

- https://academic.oup.com/bioinformatics/article/41/7/btaf375/8174968

---

# Theory and Background

## Protein-DNA Interactions

Overview of how transcription factors and chromatin-associated proteins regulate gene expression.

- https://www.nature.com/articles/nrg3682

## Histone Modifications

Review of chromatin states and histone marks.

Common marks include:

| Mark | Association |
|--------|-------------|
| H3K27ac | Active enhancers |
| H3K4me3 | Active promoters |
| H3K27me3 | Repression |
| H3K4me1 | Poised enhancers |

## CUT&RUN vs ChIP-seq

CUT&RUN generally provides:

- lower background
- higher signal-to-noise ratio
- fewer sequencing requirements
- better performance with low-input samples

compared to traditional ChIP-seq approaches.

---

# Recommended Review Articles

### Chromatin Profiling Methods

Review of approaches used to study protein-DNA interactions and chromatin states.

- https://www.nature.com/articles/s41576-018-0089-8

### CUT&RUN and CUT&Tag Technologies

Comparison of modern chromatin profiling methods.

- https://hutchdatascience.org/Choosing_Genomics_Tools/cutrun-and-cuttag.html

---

# Additional Resources

## Genome Browsers

### UCSC Genome Browser

- https://genome.ucsc.edu

### Ensembl Genome Browser

- https://www.ensembl.org

### IGV

- https://igv.org

---

## Motif Databases

### JASPAR

Open-access transcription factor motif database.

- https://jaspar.genereg.net

### HOCOMOCO

Curated vertebrate transcription factor motifs.

- https://hocomoco11.autosome.org

---

## Regulatory Databases

### ENCODE Project

- https://www.encodeproject.org

### SCREEN Regulatory Elements

- https://screen.encodeproject.org

### GeneHancer

- https://www.genecards.org/Guide/GeneHancer
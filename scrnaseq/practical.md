---
title:  Practical Session
parent: Single‑cell RNA-seq
nav_order: 2
---

# Practical Session

In this practical session, we will analyze a real single-cell RNA-seq dataset using the **Seurat** package in R.

Our goal is to follow a typical scRNA-seq workflow, starting from a count matrix and ending with the identification of distinct cell populations.

# Learning Objectives

By the end of this session, participants will be able to:

- load and inspect a single-cell RNA-seq dataset
- perform basic quality control
- normalize expression data
- identify highly variable genes
- generate low-dimensional representations of cells
- perform clustering analysis
- identify marker genes
- annotate cell populations

# Dataset

The dataset used in this exercise is the **10k Peripheral Blood Mononuclear Cells (PBMC)** dataset generated using the 10x Genomics platform.

# Workflow Overview

<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna_workflow_class.png" width="500">
</div>

<p align="center">
<em>Figure 1. Overview of the single-cell RNA-seq workflow. Source: https://hbctraining.github.io/Intro-to-scRNAseq/lessons/08_integration_cca_theory.html </em>
</p>

# About this Data

Before jumping into the data analysis, we need to have the biological context for this dataset. Here are the most relevant metadata for this dataset:

- The libraries were prepared using Chromium Next GEM Single Cell 3ʹ Reagent Kits v3.1
- The samples were sequenced on the Illumina NovaSeq 6000
- Peripheral Blood Mononuclear Cells (PBMC) samples from a healthy donor.

Since this sample is PBMCs, we will expect immune cells, such as:

- B cells
- T cells
- NK cells
- monocytes
- macrophages
- possibly megakaryocytes

Before any data analysis, this is the most important information you have to keep in mind: is this a healthy sample? If yes, what cell types were already characterized for this tissue? If not, what disease is this? Immune disease? Tumor? Diabetes? What disease markers are expected? This will save you time and tell you right away if you have a quality dataset or not.

In this particular dataset, none of the above cell types are expected to be low complexity or anticipated to have high mitochondrial content.

# Step 1 — Loading Packages

```r
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
```

# Step 2 — Set input and output paths

```r
# Set path. It's useful to do it at the beginning, so you don't have to rewrite file paths again and again
DATA_DIR <- "/course/shared/scrnaseq/data/raw"

H5_FILE <- file.path(
  DATA_DIR,
  "SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.h5"
)

OUTDIR <- "~/scrnaseq/"

dir.create(OUTDIR, showWarnings = FALSE)
```

# Step 3 — Loading Data

```r
counts <- Read10X_h5(H5_FILE)

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "PBMC10K",
  min.cells = 3,
  min.features = 200
)

# min.cells: exclude feature (genes) expressed in less than 3 cells
# min.features: exclude cells with less than 200 genes (features) expressed
```

#  About the Seurat Object
 
For bioinformatics, when we start to work with single cell, this might be the first time you've seen this type of data format. Don't panic!

The idea here is that all information related to your single cell can be accessed from one object, instead of scattered around and you running the risk of loosing it. 
 
 In a nutshell, Seurat object is an R S4 object which allows us to store single-cell data in R. 

The key slots in a Seurat object are:

- assays: This slot stores the raw and processed data in different forms. It is a list of Assay objects, each representing a specific type of data.  Examples: RNA, SCT, etc. 

- meta.data: A data.frame containing metadata associated with each cell. This can include cell type annotations, experimental conditions, or other variables related to the cells. Example columns: cell_type, batch, condition, cluster.

- reductions: A list of dimensionality reductions applied to the data. These are used for visualizations like PCA, t-SNE, or UMAP. Examples: pca, tsne, umap.

- graphs: A list of graphs (usually a nearest-neighbor graph) that are used for clustering and other analyses.  Examples: RNA_snn (a shared nearest-neighbor graph for RNA-seq data), pca_snn.

- clusters: This stores the cluster assignments for each cell after a clustering analysis (e.g., Louvain or Leiden clustering). It is typically stored in the meta.data slot but can also be stored in a separate slot.

- commands: A record of the commands used to generate the object. 

- misc: This slot is used to store arbitrary information that doesn’t fit into the other slots.

When working with single cell, the calculation can take a long time, so it's very useful to be able to do it once, save it into a new object and pull it again for other scripts.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/seurat_obj.png" width="500">
</div>

<p align="center">
<em>Figure 2. What's inside a Seurat Object. Source: https://biostatsquid.com/ </em>
</p>

```r
#Inspect Object
pbmc
```

```
An object of class Seurat 
22860 features across 10884 samples within 1 assay 
Active assay: RNA (22860 features, 0 variable features)
 1 layer present: counts
```

```r
#Inspect slots
slotNames(pbmc)
```

```
 [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"      
 [6] "neighbors"    "reductions"   "images"       "project.name" "misc"        
[11] "version"      "commands"     "tools"    
```

```r
#Inspect meta data
head(pbmc@meta.data)
```

```
                   orig.ident nCount_RNA nFeature_RNA
AAACCCAGTATATGGA-1    PBMC10K        886          343
AAACCCAGTATCGTAC-1    PBMC10K       1628          749
AAACCCAGTCGGTGAA-1    PBMC10K       6590         1867
AAACCCAGTTAGAAAC-1    PBMC10K      17318         3809
AAACCCAGTTATCTTC-1    PBMC10K       3526         1516
AAACCCAGTTGCCGAC-1    PBMC10K       6228         2110
```

```r
#Inspect data dimensions
dim(pbmc)
```

```
[1] 22860 10884
```

```r
#Inspect genes
head(rownames(pbmc))
```

```
[1] "AL627309.1" "AL627309.3" "AL627309.5" "AL627309.4" "AL669831.2" "LINC01409" 
```

```r
#Inspect cells
head(colnames(pbmc))
```

```
[1] "AAACCCAGTATATGGA-1" "AAACCCAGTATCGTAC-1" "AAACCCAGTCGGTGAA-1"
[4] "AAACCCAGTTAGAAAC-1" "AAACCCAGTTATCTTC-1" "AAACCCAGTTGCCGAC-1"
```

```r
#Inspect assays
Assays(pbmc)
```

```
[1] "RNA"
```

```r
#Inspect RNA assay
pbmc[["RNA"]]
```

```
Assay (v5) data with 22860 features for 10884 cells
First 10 features:
 AL627309.1, AL627309.3, AL627309.5, AL627309.4, AL669831.2, LINC01409,
FAM87B, LINC01128, LINC00115, FAM41C 
Layers:
 counts 
```

# Step 3 — Quality Control

```r
#Calculate Mt content
#Beware: this only works with human genome annotation!
#For other species: filter genes annotated in the mitochondrial chromosome, 
#save it to a file or vector, and use  the "features" argument to calculate de mitocondrial content correctly
pbmc[["percent.mt"]] <- PercentageFeatureSet(
  pbmc,
  pattern = "^MT-"
)
```

```r
VlnPlot(
  pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot.png" width="500">
</div>

<p align="center">
<em>Figure 3. Violin Plots. Each dot is one of the cells in our dataset. What is the overall quality of the data? </em>
</p>

</details>


```r
FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)
```


<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot01.png" width="500">
</div>

<p align="center">
<em>Figure 4. Scatter plot of RNA count per cell vs MT content. What pattern do you see? </em>
</p>

</details>


```r
FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot02.png" width="500">
</div>

<p align="center">
<em>Figure 5. Scatter plot of RNA count and RNA feature What does this tell us about the quality? </em>
</p>

</details>


# Step 4 — Filtering Cells

```r
#Filter Based on QC
pbmc <- subset(
  pbmc,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 5000 &
    percent.mt < 10
)

#Verify dims
pbmc

dim(pbmc)
```

How many cells did we eliminate? Why did we choose these thresholds?

# Step 5— Data Normalization

```r
#To remove technical variations, sequencing depth and variance amongst our cells, we log normalize our counts 
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)
```

# Step 6 — Identifying Highly Variable Genes

```r
#Find Highly Variable Genes
pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

#Visualize
plot1 = VariableFeaturePlot(pbmc)

top10 <- head(VariableFeatures(pbmc), 10)

top10

LabelPoints(plot = plot1, points = top10, repel = TRUE)

```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot03.png" width="500">
</div>

<p align="center">
<em>Figure 6. Scatter Plot of Average Expression versus Standardized Variance. What could you tell about this scattering? What type of genes are the most variable? </em>
</p>

</details>


  <details>
<summary>Click to reveal table</summary>

<table>
<thead>
<tr>
<th>Gene</th>
<th>Biological Function</th>
<th>Likely Cell Type</th>
</tr>
</thead>
<tbody>
<tr><td><b>IGLC1</b></td><td>Immunoglobulin lambda constant 1</td><td>B cells / Plasma cells</td></tr>
<tr><td><b>IGLC3</b></td><td>Immunoglobulin lambda constant 3</td><td>B cells / Plasma cells</td></tr>
<tr><td><b>CXCL10</b></td><td>Interferon-induced inflammatory chemokine</td><td>Activated immune cells / Interferon-responsive cells</td></tr>
<tr><td><b>PTGDS</b></td><td>Prostaglandin D2 synthase</td><td>Dendritic cells / pDCs</td></tr>
<tr><td><b>GZMB</b></td><td>Granzyme B</td><td>Cytotoxic NK cells</td></tr>
<tr><td><b>JCHAIN</b></td><td>Joining chain of IgA and IgM antibodies</td><td>Plasmablasts / Plasma cells</td></tr>
<tr><td><b>CDKN1C</b></td><td>Cell cycle inhibitor</td><td>Monocyte/DC subsets</td></tr>
<tr><td><b>GNLY</b></td><td>Granulysin</td><td>NK cells</td></tr>
<tr><td><b>MZB1</b></td><td>Antibody secretion chaperone</td><td>Plasmablasts</td></tr>
<tr><td><b>LILRA4</b></td><td>Leukocyte immunoglobulin-like receptor A4</td><td>Plasmacytoid dendritic cells (pDCs)</td></tr>
</tbody>
</table>

</details>


Why identify highly variable genes? For huge datasets, this step can cut off a lot of memory usage in downstream calculations. This is also a good time to check the quality of the data. Ask yourself: What kind of genes appear here? Do they have biological relevance? Or is it ribosomal or mitochondrial genes? If that's the case, this is a good time to go back and review your QC: you might have been afraid of "losing" cells, but kept all the bad apples.

# Step 7 — Principal Component Analysis (PCA)

Principal Components Analysis is a dimensionality-reduction method employed to stratify data by their variance. In other words, points that are more closely together on a PCA are more similar to each other, while points that are more dissimilar to each other will be further apart. 

This is very useful in any transcriptomic data, where you have the expression value for every gene in multiple or thousands of samples/cells. We are talking about 25,000 or more genes expressed. It would be very laborious to try and plot all possible combinations of two genes to extract any information at all of how similar two samples or cells are. Luckily, we don't have to suffer what others have suffered before and resolved already: we can use PCA to extract this information.

```r
all.genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all.genes
)

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(pbmc)
)

#Visualize
DimPlot(
  pbmc,
  reduction = "pca"
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot10.png" width="500">
</div>

<p align="center">
<em>Figure 7. Dimesion plot of the first 2 PCs </em>
</p>

</details>

```r
ElbowPlot(pbmc)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot04.png" width="500">
</div>

<p align="center">
<em>Figure 8. Elbow Plot: how much each PC represents the data (in percentage) </em>
</p>

</details>

# Step 8 — UMAP Visualization

As you can see from our previous plot, each PC can only explain so much about the data at a time. A way to try and use as much information as possible to cluster our cells is to employ an Uniform Manifold Approximation and Projection (UMAP). While PCA will determine all PCs, we can only plot two at a time. In contrast, UMAP will take the information from any number of top PCs to arrange the cells in this multidimensional space. 

```r
#Calculate UMAP
pbmc <- RunUMAP(
  pbmc,
  dims = 1:20
)

#Visualize
DimPlot(
  pbmc,
  reduction = "umap"
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot05.png" width="500">
</div>

<p align="center">
<em>Figure 8. UMAP calculated using 20 PCs </em>
</p>

</details>

# Step 9 — Cell Clustering

Now that we have our UMAP, we will use Seurat's graph-based clustering approach using a K-nearest neighbor to from our clusters based on distance between dots, which here is a proxy for similarity of the transcriptome of the cells.

<div align="center">
<img src="/genomics-course-usp2026/assets/images/clustering_example.png" width="500">
</div>

<p align="center">
<em>Figure 9. Schematic representation of K-nearest neighbor–based clustering for single-cell RNA-seq data. source: https://hbctraining.github.io/Intro-to-scRNAseq/lessons/10_clustering_cells_SCT.html </em>
</p>

```r
#Calculate Clusters
pbmc <- FindNeighbors(
  pbmc,
  dims = 1:20
)

pbmc <- FindClusters(
  pbmc,
  resolution = 0.5
)

DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot06.png" width="500">
</div>

<p align="center">
<em>Figure 10. UMAP with calculated  clusters </em>
</p>

</details>

# Saving

Before moving on, this is a great time to save your processed object:

```r
#Save Clustered Data
saveRDS(
  pbmc,
  file.path(
    OUTDIR,
    "pbmc_clustered.rds"
  )
)
```

This object now can be read in another scripts using the follow command:

```r
#Open RDS object

DATA_DIR <- "~/scrnaseq/"

RDS_FILE <- file.path(
  DATA_DIR,
  "pbmc_clustered.rds"
)

pbmc_processed  = readRDS(
  file.path(
    OUTDIR,
    "pbmc_clustered.rds"
  )
)
```

Remember you can do this at any step! I highly recommend it after filtering and after UMAP+clustering. Save your future self the trouble of having to rerun a complete qc/filtering/normalization/pca/umap anytime you want to do a simple feature plot. This will also help your results to be more consistent and reproducible.

# Step 10 — Identifying Marker Genes

```r
markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

#Select top10 markers per cluster
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

#Inspect
head(top10)

#Save to CSV
write.csv(
  markers,
  file.path(
    OUTDIR,
    "cluster_markers.csv"
  ),
  row.names = FALSE
)
```

Here I also recommend saving the csv file for easiness of analysis downstream. This will make your following scripts run faster. The csv file can be easily manipulated outside R.

# Step 11 — Cell Type Annotation

## Explore some known marker

```r
FeaturePlot(pbmc, features = "MS4A1")
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot07.png" width="500">
</div>

<p align="center">
<em>Figure 11. UMAP showing expression of Gene MS4A1 </em>
</p>

</details>

```r
FeaturePlot(pbmc, features = "CD3D")
FeaturePlot(pbmc, features = "NKG7")
FeaturePlot(pbmc, features = "LYZ")
```

## Explore top markers per cluster

```r
# Let's use our top markers to infer cell identity
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  ungroup()

DoHeatmap(
  pbmc,
  features = unique(top5$gene),
  size = 3
) +
  NoLegend()
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot08.png" width="700">
</div>

<p align="center">
<em>Figure 12. Heatmap of top 5 markers across clusters. Do you recognize any of them from our previous images? </em>
</p>

</details>

## Explore known markers all at the same time

```r
# And use the canonical markers
canonical_markers <- list(
  "T cells" = c("CD3D", "IL7R"),
  "NK cells" = c("NKG7", "GNLY"),
  "B cells" = c("MS4A1", "CD79A"),
  "Monocytes" = c("LYZ"),
  "Dendritic cells" = c("FCER1A", "LILRA4"),
  "Plasma cells" = c("JCHAIN")
)

DotPlot(
  pbmc,
  features = canonical_markers
) +
  RotatedAxis()
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot09.png" width="700">
</div>

<p align="center">
<em>Figure 13. Dotplot of selected markers across clusters. Dot size represents percentage of cells from a given cluster that express the gene </em>
</p>

</details>

Given the markers, we can make a best guest on the cell types:

```r
0,2,8,11 → Monocytes, or maybe more than one type ofmonocytes?
6,7,14 → B cells
5,9,10 → NK cells
12 → Dendritic cells / pDC-like
1,3,4 → T cells
```

To annotate our cells, we add a column to the meta.data

```r
cluster_annotations <- c(
  "0"  = "Monocytes",
  "1"  = "T cells",
  "2"  = "Monocytes",
  "3"  = "T cells",
  "4"  = "T cells",
  "5"  = "NK cells",
  "6"  = "B cells",
  "7"  = "B cells",
  "8"  = "Monocytes",
  "9"  = "NK cells",
  "10" = "NK cells",
  "11" = "Monocytes",
  "12" = "Dendritic cells",
  "13" = "Monocytes",
  "14" = "B cells"
)

pbmc_$celltype_manual <- unname(
  cluster_annotations[
    as.character(Idents(pbmc))
  ]
)

```

And now we inspect the annotated UMAP

```r
DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE
)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/scrna/Rplot11.png" width="500">
</div>

<p align="center">
<em>Figure 14. UMAP of Annotated clusters </em>
</p>

</details>

Would you be satisfied with this annotation? How could you improve it?

Finally, save our annotated object:

```r
saveRDS(
  pbmc,
  file.path(
    OUTDIR,
    "pbmc_annotated.rds"
  )
)
```

# The End

I know you would like to keep analyzing this single cell data to it's fullest and discover all the subclasses of the known cell types and publish your very own Nature paper. But fear not! Though this class has come to an end, I'll leave resources in the  [References](./references)  tab to go even deeper and wilder in your analysis:
- what happens when I have a healthy and a disease sample?
- what about control versus treatment?
- what if I have different stages of the same tissue or disease?

There's also many public datasets available and you can use them to explore your own hypothesis and integrate them into your project.

Go on and become the single cell specialist that your lab desperately needs!




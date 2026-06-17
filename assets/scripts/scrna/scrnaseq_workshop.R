#Load Packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)

#Set Paths
DATA_DIR <- "~/course/shared/scrnaseq/data/raw"

H5_FILE <- file.path(
  DATA_DIR,
  "SC3_v3_NextGem_SI_PBMC_10K_filtered_feature_bc_matrix.h5"
)

OUTDIR <- "~/scrnaseq"

dir.create(OUTDIR, showWarnings = FALSE)

#Load Data
counts <- Read10X_h5(H5_FILE)

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "PBMC10K",
  min.cells = 3,
  min.features = 200
)

#Inspect Data
pbmc

slotNames(pbmc)

head(pbmc@meta.data)

dim(pbmc)

head(rownames(pbmc))

head(colnames(pbmc))

Assays(pbmc)

pbmc[["RNA"]]

#Calculate Mt content
#Beware: this only works with human genome annotation!
pbmc[["percent.mt"]] <- PercentageFeatureSet(
  pbmc,
  pattern = "^MT-"
)

#Visualize QC
VlnPlot(
  pbmc,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

FeatureScatter(
  pbmc,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

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

#Save Filtered Object
saveRDS(
  pbmc,
  file.path(
    OUTDIR,
    "pbmc_qc.rds"
  )
)

#Normalize
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

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

#Calculate PCA
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

ElbowPlot(pbmc)

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

#Save Clustered Data
saveRDS(
  pbmc,
  file.path(
    OUTDIR,
    "pbmc_clustered.rds"
  )
)

#Find Markers
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

#Inspect for Gene Markers
#B Cells
FeaturePlot(pbmc,
            features = "MS4A1")

#T Cells
FeaturePlot(pbmc,
            features = "CD3D")

#NK Cells
FeaturePlot(pbmc,
            features = "NKG7")
FeaturePlot(pbmc, features="GNLY")
FeaturePlot(pbmc, features="GZMB")

#Monocytes
FeaturePlot(pbmc,
            features = "LYZ")

#Dendeitic cells
FeaturePlot(pbmc,
            features = "FCER1A")

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

pbmc$celltype_manual <- cluster_annotations[
  as.character(Idents(pbmc))
]

DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE
)

#save object
saveRDS(
  pbmc,
  file.path(
    OUTDIR,
    "pbmc_annotated.rds"
  )
)
---
title: "scRNA-SEQ-course-2019"
author: "David John & Wesley Abplanalp"
date: "12 6 2019"
---


# 1.) Install required Software
## 1.2.) Install R
1. Windows: <https://cran.rstudio.com/bin/windows/base/>
2. Mac:<https://cran.rstudio.com/bin/macosx/>
3. Linux: Please run the following commands in the terminal
```{shell}
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

sudo apt update

sudo apt install r-base r-base-core r-recommended
```

## 1.2.) Install Rstudio
To install RStudio please follow the link and install the right RStudio for your operating system.

1. Windows: <https://download1.rstudio.org/desktop/windows/RStudio-1.2.1335.exe>
2. Mac:<https://download1.rstudio.org/desktop/macos/RStudio-1.2.1335.dmg>
3. Linux: <https://www.rstudio.com/products/rstudio/download/#download>

## 1.3.) Install required R-packages
For or training we will need the following R packages, which can be installed by the following R commands
```{r}
# Enter commands in R (or R studio, if installed)
# Install the devtools package from Hadley Wickham
install.packages('devtools')
# Replace '2.3.0' with your desired version
devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))

install.packages("tidyr")
install.packages("dplyr")
install.packages("ggpubr")
install.packages("scales")
install.packages("ggplot2")
install.packages("stringr")

```

## 1.4.) Test the required R-packages
To test all installed packages, run each command individually. 
If the installed package are loaded witout any error message, the package was installed sucessfully. 
```{r}
library(tidyr)
library(dplyr)
library(ggpubr)
library(scales)
library(stringr)
library(Seurat)
library(ggplot2)
```

#2.) Download The required Datasets
During this course we will analyse 3 datasets. Two are punlished healthy PBMC's from 10x and one is our data, from a heart failure patient of the university clinic Frankfurt.
1.) 2700 PBMCs (v1 Chemistry, healthy) 
*(Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k>)
2.) 5000 PBMCs (V3 Chemistry, healthy) 
*(Details: <https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_v3>)
3.) 1000 CD31 positive cells 
*(V2 Chemistry, heart failure)


##2.1.) Load the Datasets
For this course, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. 
There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. The raw data can be found here.

Detailed informations about this dataset can be found here:



##### TODO:
- [ ] Download the datasets, save it to your PC and unzip the file:
<http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>


##2.2) Load the RScript
-[ ] Download the Rscript (<https://github.com/djhn75/scRNA-SEQ-ZMM/blob/master/scRNA-seq-course-2019.Rmd>) 
-[ ] Open the file in Rstudio and continue at step 3



#3.) Analyze 2700 PBMCs)

##3.1.) Import to Seurat
```{r}
library(dplyr)
library(Seurat)

# Load the PBMC dataset
tmp.object <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc3k_v1 <- CreateSeuratObject(counts = tmp.object, project = "pbmc3k_v1", min.cells = 3, min.features = 200)
```

##3.1.) Pre-processing The Raw data

###3.1.1.) QC and selecting cells for further analysis
The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.
QC and selecting cells for further analysis

Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics commonly used by the community include
  
  *The number of unique genes detected in each cell.
    *Low-quality cells or empty droplets will often have very few genes
    *Cell doublets or multiplets may exhibit an aberrantly high gene count
  *Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
  *The percentage of reads that map to the mitochondrial genome
    *Low-quality / dying cells often exhibit extensive mitochondrial contamination
    *We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
    *We use the set of all genes starting with MT- as a set of mitochondrial genes

```{r}
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc3k_v1[["percent.mt"]] <- PercentageFeatureSet(pbmc3k_v1, pattern = "^MT-")
```
TODO: __Where are QC metrics stored in Seurat?__

In the example below, we visualize QC metrics, and use these to filter cells.

  *We filter cells that have unique feature counts over 2,500 or less than 200
  *We filter cells that have >5% mitochondrial counts
```{r}
# Visualize QC metrics as a violin plot
VlnPlot(pbmc3k_v1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc3k_v1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc3k_v1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
pbmc3k_v1 <- subset(pbmc3k_v1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```

###3.1.2.) Normalizing the data 
After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in `pbmc[["RNA"]]@data`.

```{r}
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```


###3.1.3.) Identification of highly variable features (feature selection)
We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.

Our procedure in Seurat3 is described in detail here, and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the `FindVariableFeatures` function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.
```{r}
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
```


###3.1.3.) Scaling the data
Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData function:

  *Shifts the expression of each gene, so that the mean expression across cells is 0
  *Scales the expression of each gene, so that the variance across cells is 1
    *This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
  *The results of this are stored in `pbmc[["RNA"]]@scale.data`
```{r}
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
```


##3.2.) Clustering 
###3.2.1.) Perform linear dimensional reduction
Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.

Seurat provides several useful ways of visualizing both cells and features that define the PCA, including `VizDimReduction`, `DimPlot`, and `DimHeatmap`

```{r}
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```

An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one (`ElbowPlot` function). In this example, we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs.
```{r}
ElbowPlot(pbmc)

```
Identifying the true dimensionality of a dataset – can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-12 as a cutoff.

We chose 10 here, but encourage users to consider the following:

    *Dendritic cell and NK aficionados may recognize that genes strongly associated with PCs 12 and 13 define rare immune subsets (i.e. MZB1 is a marker for plasmacytoid DCs). However, these groups are so rare, they are difficult to distinguish from background noise for a dataset of this size without prior knowledge.
    *We encourage users to repeat downstream analyses with a different number of PCs (10, 15, or even 50!). As you will observe, the results often do not differ dramatically.
    *We advise users to err on the higher side when choosing this parameter. For example, performing downstream analyses with only 5 PCs does signifcanltly and adversely affect results.


###3.2.2.) Cluster the cells
Seurat v3 applies a graph-based clustering approach, building upon initial strategies in (Macosko et al). Importantly, the distance metric which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [SNN-Cliq, Xu and Su, Bioinformatics, 2015] and CyTOF data [PhenoGraph, Levine et al., Cell, 2015]. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.

As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the FindNeighbors function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).

To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [SLM, Blondel et al., Journal of Statistical Mechanics], to iteratively group cells together, with the goal of optimizing the standard modularity function. The FindClusters function implements this procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the Idents function.
```{r}
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

```{r}
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```



###3.2.3.) Run non-linear dimensional reduction (UMAP/tSNE)
Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.
```{r}
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

```{r}
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
```

#####Save
You can save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above, or easily shared with collaborators.
```{r}
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
```



###3.2.4.) Finding differentially expressed features (cluster biomarkers)
Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed features will likely still rise to the top.

```{r}
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
```

```{r}
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
```

```{r}
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
```


###3.2.5 Plot DEG's
We include several tools for visualizing marker expression. VlnPlot (shows expression probability distributions across clusters), and FeaturePlot (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot, CellScatter, and DotPlot as additional methods to view your dataset.

```{r}
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
```

```{r}
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
```

```{r}
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
    "CD8A"))
```

DoHeatmap generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

```{r}
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```


###3.2.1.) Assigning cell type identity to clusters
Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:
```{r}
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
```

#####Save
```{r}
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")
```


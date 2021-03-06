library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)
require(scales)
library(ggpubr)
library(tidyr)


#' 
#' @author David John
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCellsByQuantile <- function(seuratObject, lowQuantile=0.1 , highQuantile=0.95, maxMito=0.1){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  cat("FilterByQuantile\n")
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  lowQuantile<<-lowQuantile
  highQuantile<<-highQuantile
  maxMito<<-maxMito
  sample<-unique(seuratObject$sample)
  Quality <- data.frame(UMI=seuratObject$nCount_RNA, nGene=seuratObject$nFeature_RNA, label = factor(seuratObject$sample), percent.mito=seuratObject$percent.mito)
  
  Quantile.low.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = lowQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.high.UMI <- Quality %>% group_by(label) %>%
    summarise(UMI = list(enframe(quantile(UMI,probs = highQuantile)))) %>%
    unnest(cols = c(UMI))
  
  Quantile.low.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = lowQuantile)))) %>%
    unnest(cols = c(nGene))
  
  Quantile.high.Gene <- Quality %>% group_by(label) %>%
    summarise(nGene = list(enframe(quantile(nGene,probs = highQuantile)))) %>%
    unnest(cols = c(nGene))
  
  
  df<-seuratObject@meta.data
  
  gg1<- ggplot(Quality, aes(x="nUMI", y=UMI)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.UMI$value, label=Quantile.high.UMI$value , vjust = -1)) + 
    geom_hline(yintercept = Quantile.low.UMI$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.UMI$value, label=Quantile.low.UMI$value , vjust = -1))
  
  gg2<- ggplot(Quality, aes(x="nFeature_RNA", y=nGene)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = Quantile.high.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.high.Gene$value, label=Quantile.high.Gene$value , vjust = -1)) +   geom_hline(yintercept = Quantile.low.Gene$value, color="red", linetype="dashed") + geom_text(aes(0.9,Quantile.low.Gene$value, label=Quantile.low.Gene$value , vjust = -1))
  
  
  gg3<- ggplot(Quality, aes(x=" % Mt Content", y=percent.mito)) + geom_violin(scale = "width") + 
    theme(axis.title.x = element_blank(),axis.ticks.x = element_blank(), legend.position = "none", axis.text.x = element_text(size=12, face = "bold"), axis.title.y = element_blank(), axis.text.y = element_text(size=10)) + 
    geom_hline(yintercept = maxMito, color="red", linetype="dashed") + geom_text(aes(0.9,maxMito, label=maxMito , vjust = -1))
  
  gg<-ggarrange(gg1,gg2,gg3, ncol = 3)
  
  library(ggpubr)  
  
  gg<-annotate_figure(gg, fig.lab = sample, fig.lab.pos = "top", fig.lab.size = 15, fig.lab.face = "bold")
  
  seuratObject<- subset(x= seuratObject, subset = nCount_RNA < Quantile.high.UMI$value & nCount_RNA > Quantile.low.UMI$value & 
                          nFeature_RNA < Quantile.high.Gene$value & nFeature_RNA > Quantile.low.Gene$value & percent.mito < maxMito)
  
  
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  cat("Filtered ",diff, "from" , sizeBefore, " cells\n", "(minFeatures=",Quantile.low.Gene$value, "; maxFeatures=", Quantile.high.Gene$value, "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(maxMito)
  return(list(seuratObject, gg))
}





#' Import Single cell sequencing experiments into Seurat3and perform normalisation and scale Data 
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE, FilterCells=TRUE, FilterByAbsoluteValues=FALSE) {
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  if (length(mito.features)<10) {
    mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratObject), value = TRUE)
  }
  if (length(mito.features)<10) {
    stop("Error: Could not find MT genes")
  }
  
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts'))
  seuratObject$percent.mito <- percent.mito
  
  #write QC to file
  svg(paste0(pathway,"QC_preFiltered.svg"))
  gg_preFiltering<-VlnPlot(object = seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0,)
  print(gg_preFiltering)
  dev.off()
  cat("start Filtering")
  if (FilterCells==TRUE) {
    if (FilterByAbsoluteValues==TRUE) {
      seuratObject<-FilterDeadCells(seuratObject = seuratObject)
    }
    else {
      tmp<-FilterDeadCellsByQuantile(seuratObject = seuratObject, lowQuantile = 0.1, highQuantile = 0.95)
      seuratObject<-tmp[[1]]
      svg(paste0(pathway,"QC_QuantileFiltering.svg"))
      print(tmp[[2]])
      dev.off()
      gg_preFiltering<-tmp[[2]]
      
    }
    
  }
  if (performNormalisation==TRUE) {
    seuratObject<-NormalizeData(object = seuratObject,verbose = FALSE)
  }
  if(performVariableGeneDetection==TRUE){
    seuratObject<-FindVariableFeatures(object = seuratObject, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  if (performScaling==TRUE) {
    seuratObject<-ScaleData(object = seuratObject)
  }
  cat("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
  return(list(seuratObject, gg_preFiltering))
}


#' 
#' @author David John
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCells <- function(seuratObject, minFeatures=200, maxFeatures=3000, maxMito=0.1){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  sizeBefore<-length(seuratObject@meta.data$orig.ident)
  
  #For some unknown reasons these variables need to be global for the subset function, otherwise there is an eval() unknown variable error 
  minFeatures<<-minFeatures
  maxFeatures<<-maxFeatures
  maxMito<<-maxMito
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > minFeatures & nFeature_RNA < maxFeatures & percent.mito < maxMito)
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  cat("Filtered ",diff, "from" , sizeBefore, " cells\n", "(minFeatures=",minFeatures, "; maxFeatures=", maxFeatures, "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
  rm(minFeatures,maxFeatures,maxMito)
  return(seuratObject)
}

#' Import and combine several Single cell sequencing experiments into Seurat
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineSeuratObjects <- function(pathways,ids, performNormalisation = FALSE, performVariableGeneDetection=FALSE){
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  for (i in 1:length(pathways)) {
    if (i<2) { 
      seuratObject1<-Importer(pathways[i],ids[i], TenX = TRUE, performNormalisation = performNormalisation, performVariableGeneDetection = performVariableGeneDetection)
      next
    }
    seuratObject2<-Importer(pathways[i],ids[i])
    seuratObject1 <- merge(x = seuratObject1, y = seuratObject2) 
  }
  cat("Merged Seurat object contains ", length(seuratObject1@meta.data$orig.ident)," cells\n")
  return(seuratObject1)
}



#library(scImpute)
library(Seurat)
#' Read 10X results and impute the dataset with scImpute
#' @author David John
#' @param pathways Pathway to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids String that are assigned to the output matrix
imputeData <- function(pathways,ids, cluster=12, ncores=20, drop_thre=0.5){
  for (i in 1:length(pathways)) {
    cat("Start imputing Sample ", pathways[i] )
    path.Matrix<-paste(pathways[i],"Matrix.csv",sep="/")
    path.Imputed.Matrix <- paste(pathways[i], "scImpute", ids[i], sep="/")
    Ten_X <- Importer(pathways[i], id = ids[i], TenX = TRUE, performNormalisation = FALSE)
    Ten_X <- FilterDeadCells(seuratObject = Ten_X)
    cat("Write temporary Martix to ", path.Matrix)
    write.csv(as.data.frame(as.matrix(Ten_X@assays$RNA@data)), file = path.Matrix)
    cat("Start imputation for ", path.Matrix)
    #scimpute(count_path = path.Matrix, out_dir = path.Imputed.Matrix, Kcluster = cluster, ncores=ncores, drop_thre = drop_thre)
    cat("Wrote imputed Martix to ", path.Imputed.Matrix)
  }
}



#' Import and combine imputet single cell Matrices into Seurat with the new method VST of Seurat3
#' @author David John
#' @param pathways A vector of pathways to scimputed count marices 
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineScImputedSeuratObjectsVST <- function(pathways,ids){
  
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  hvg<-c() #to save high variable genes
  seurat.objects<-list() #list to dave Seurat objects
  #loop through directories of pathways
  for (i in 1:length(pathways)) {
    matrixPath<-paste(pathways[i],"scImpute",paste(ids[i],"scimpute_count.csv",sep=""),sep = "/")
    seuratObject<-Importer(matrixPath,id = ids[i], TenX = FALSE)
    seurat.objects<-c(seurat.objects,seuratObject)
  }
  
  cat("start Integration")
  seurat.objects.anchors <- FindIntegrationAnchors(object.list = seurat.objects)
  seuratObject <- IntegrateData(anchorset = seurat.objects.anchors)
  seuratObject <- ScaleData(object = seuratObject, verbose = FALSE)
  cat("Merged Seurat object contains ", length(seuratObject@meta.data$orig.ident)," cells\n")
  
  return(seuratObject)
}



#' Drawa FeatureHeatmap
#' @author Lukas Tombor and David John
#' @param object A Seurat3 object
#' @param features Genes which are plotted
#' @param group.by Split/Group TSNES by this parameter
#' @param cols Colors of visible cells
FeatureHeatmap <- function(object, features, group.by=NULL, cols=c("skyblue","red4"), assay="RNA") {
  require(reshape2);require(ggplot2)
  DefaultAssay(object)<-"RNA"
  if (!group.by %in% colnames(object@meta.data)) {
    stop("Grouping parameter was not found in the meta data")
  }
  
  #test if Genes can be found
  for(i in 1:length(features)){
    if (!features[i] %in% rownames(object@assays$RNA@data)) {
      stop("Gene was not found in the Expression Matrix")
    }
  }
  
  
  A <- data.frame(object@meta.data)
  X <- Embeddings(object = object, reduction = "tsne")
  coord = NULL
  for(i in rownames(A)){
    coord <- rbind(coord, c(X[i,], i))
  }
  nclusters<-length(table(A[,group.by]))
  A <- data.frame(A, coord)
  A$tSNE_1 <- as.numeric(levels(A$tSNE_1)[A$tSNE_1])
  A$tSNE_2 <- as.numeric(levels(A$tSNE_2)[A$tSNE_2])
  A$seurat_clusters <- factor(A$seurat_clusters, levels = 0:(nclusters-1))
  #A$sample_dummy <- as.integer(as.factor(A$sample))
  
  for(i in 1:length(features)){
    a.colnames<-colnames(A)
    a.colnames<-c(a.colnames,paste0(features[i],".exprs"))
    A <- data.frame(A, x=GetAssayData(object, assay = assay)[features[i], ])
    colnames(A)<-a.colnames
  }
  
  A.rep<-as.data.frame(lapply(A, rep, nclusters))
  f3 <- function(variables) {
    return(names(table(A.rep[,group.by])))
  }
  
  A.rep$Cluster_fate <- ave(rep(NA,nrow(A.rep)), A.rep[,group.by], FUN = f3)
  
  #add new column with NA for unfitting cluster and expression value for the correct cluster
  for (f in features) {
    A.rep[,f]<-ifelse(A.rep[,group.by]==A.rep$Cluster_fate,A.rep[,paste0(f,".exprs")],NA)
  }
  
  A.melt<-melt(A.rep, measure.vars = colnames(A.rep[,c((ncol(A.rep)-1):ncol(A.rep))]))
  
  ggplot(A.melt, aes(x=tSNE_1, y= tSNE_2, color = value))+
    geom_point(size = 0.2)+
    scale_color_continuous(low = cols[1], high = cols[2], na.value = "grey93", name="Scaled Expression")+
    facet_grid(c("variable","Cluster_fate"), scales = "free", drop = FALSE)+
    labs(y="Genes", x= group.by, title = "")+
    theme_bw()+
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "top", strip.text = element_text(face = "italic"))
  
}


PlotClusterPrecentageBarplot <- function(x, faceWrapeBy, fileName=NULL, width=15, height=8){
  df<-data.frame(x1=x,x2=faceWrapeBy)
  
  df <- df %>% count(x1,x2,.drop = FALSE) %>% group_by(x1) %>%
    mutate(freq = n /sum(n)) %>% complete(x2,fill = list(n=0,freq=0))
  if (!is.null(fileName)) {
    svg(filename = fileName,width = width, height = height)
    
  }
  ggplot(df, aes(x=x1, y= freq, fill= x1))+
    geom_col(width = 0.9, color = "black")+
    facet_wrap(~x2, nrow = 2, scales = "free")+
    scale_y_continuous(name = "Percent per timepoint", labels = scales::percent_format())+
    theme(panel.background = element_blank(), axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust= 1, size = 8))
  
  
  
}

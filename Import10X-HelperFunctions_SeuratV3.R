library(Seurat)
library(dplyr)
library(tibble)
library(reshape2)
require(scales)
library(ggpubr)
library(tidyr)
#library(EnhancedVolcano)



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





#' Import Single cell sequencing experiments into Seurat3and perform normalisation and scale Data and do a summary of mapping stats, optional perform Doubletfinder
#' @author David John & Mariano Ruz Jurado
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performScaling = FALSE,performVariableGeneDetection=TRUE, FilterCells=TRUE, FilterByAbsoluteValues=FALSE,...) {
  require(Seurat)
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  #catch optional parameters 
  optionalParameters <- list(...)
  
  seuratObject =CreateSeuratObject(counts = Matrix, project = id, min.cells = 5)
  seuratObject$sample <- id
  tmp<-unlist(strsplit(id,split = "-"))
  seuratObject$condition <- paste0(tmp[1:length(tmp)-1],collapse = "-")
  
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  if (length(mito.features)<10) {
    mito.features <- grep(pattern = "^mt-", x = rownames(x = seuratObject), value = TRUE)
  }
  if (length(mito.features)<1) {
    warning("Error: Could not find MT genes")
    
  }
  
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts'))
  seuratObject$percent.mito <- percent.mito
  
  #write QC to file
  p1<-VlnPlot(object = seuratObject, features = c("nFeature_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p2<-VlnPlot(object = seuratObject, features = c("nCount_RNA"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  p3<-VlnPlot(object = seuratObject, features = c("percent.mito"), ncol = 1, pt.size = 0) + theme(legend.position = "None", axis.title.x = element_blank(), axis.text.x = element_blank())
  gg_preFiltering <- ggarrange(p1,p2,p3, nrow = 1)
  gg_preFiltering <- annotate_figure(gg_preFiltering, top = text_grob(id,face="bold",color = "darkred",size=18,hjust = 0.2))
  ggsave(filename = paste0(pathway,"QC_preFiltered.svg"),device = "svg", width = 10,height = 10)
  
  if (FilterCells==TRUE) {
    message("start Filtering")
    if (FilterByAbsoluteValues==TRUE) {
      if (is.null(optionalParameters$minFeatures)) {
        stop("Please define 'minFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxFeatures)) {
        stop("Please define 'maxFeatures' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$minCounts)) {
        stop("Please define 'minCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxCounts)) {
        stop("Please define 'maxCounts' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      if (is.null(optionalParameters$maxMito)) {
        stop("Please define 'maxMito' while filtering for absolute values (FilterByAbsoluteValues==TRUE)")
      }
      message("Running FilterDeadCells")
      seuratObject<-FilterDeadCells(seuratObject = seuratObject,
                                    minFeatures = optionalParameters$minFeatures,
                                    maxFeatures = optionalParameters$maxFeatures,
                                    minCounts = optionalParameters$minCounts,
                                    maxCounts = optionalParameters$maxCounts,
                                    maxMito = optionalParameters$maxMito)
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
  message("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
  
  
  return(list(seuratObject, gg_preFiltering))
}

#' @author Mariano Ruz Jurado
#' @param pathway A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param id Vector of strings that are assigned to the concordant cells
#' @return a list with summary about mapping results 
SummarizeMapping <- function(pathway,id){
## mapping stats, i hope everyone uses cellranger and starsolo directories for these analysis, else no summary
if (file.exists(paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv"))) {
  metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "outs"))[1],"outs/metrics_summary.csv")
  #define your own numeric class
  setClass('myNum')
  #define conversion
  setAs("character", "myNum",
        function(from) as.numeric(gsub(",","", gsub("%","",from))))
  #read data with custom colClasses
  metrics_summ <- read.csv(metrics_summ.path,
                       header = T,   
                       stringsAsFactors=FALSE,
                       colClasses=c("myNum"))
  
  typeof(metrics_summ$Fraction.Reads.in.Cells)
  
  metrics_col <- as.data.frame(colnames(metrics_summ))
  rownames(metrics_col) <- metrics_col[,1]
  metrics_col[,1] <- as.character(as.vector(metrics_summ[1,]))
  metrics_summ <- metrics_col
  
  # warnings CELLRANGER
  if (metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),] < 70) {
    warning(paste0(id,": Reads mapped confidently to genome only ", metrics_summ[grep(pattern = "Confidently.to.Genome",rownames(metrics_summ)),]))
  }
  if (metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),] < 30) {
    warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Confidently.to.Transcriptome",rownames(metrics_summ)),]))
  }
  if (paste(unlist(strsplit(metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),],split=",")),collapse = "") < 1000) {
    warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number.of.Cells",rownames(metrics_summ)),]), " ,maybe the 0s were cut because of CR way of displaying numbers,  if unsure check CR web_summary")
  }
  if (as.numeric(paste(unlist(strsplit(metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),],split=",")),collapse = "")) < 300) {
    warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median.Genes.per.Cell",rownames(metrics_summ)),])," ,maybe the 0s were cut because of CR way of displaying numbers, if unsure check CR web_summary")
  }
} else {
  metrics_summ.path <- paste0(unlist(strsplit(pathway,split = "Gene"))[1],"Gene/Summary.csv")
  metrics_summ <- read.delim2(metrics_summ.path, header = F, sep = ",")
  rownames(metrics_summ) <- metrics_summ[,1]
  metrics_summ[,1] <- NULL
  
  # warnings STAR
  if (metrics_summ[7,] < 0.70) { # mapped to genome, no grep since same name as other row 
    warning(paste0(id,": Reads mapped confidently to genome only ",metrics_summ[7,]))
  }
  if (metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),] < 0.30) {
    warning(paste0(id,": Reads mapped confidently to transcriptome only ", metrics_summ[grep(pattern = "Transcriptome: Unique Genes",rownames(metrics_summ)),]))
  }
  if (metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),] < 1000) {
    warning(paste0(id,": Estimated Number of Cells only ", metrics_summ[grep(pattern = "Number of Cells",rownames(metrics_summ)),]))
  }
  if (metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),] < 300) {
    warning(paste0(id,": Median Genes per Cell only ", metrics_summ[grep(pattern = "Median Genes per Cell",rownames(metrics_summ)),]))
  }
} 
return(metrics_summ)
}

#' 
#' @author David John & Mariano Ruz Jurado
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCells <- function(seuratObject, minFeatures=300, maxFeatures=6000,minCounts=500,maxCounts=15000, maxMito=0.05){
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
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > minFeatures & nFeature_RNA < maxFeatures & nCount_RNA > minCounts & nCount_RNA < maxCounts & percent.mito < maxMito)
  
  diff<-  sizeBefore -length(seuratObject@meta.data$orig.ident)
  percent <- round(diff/sizeBefore*100,digits = 2)
  cat("Filtered ",diff, "from" , sizeBefore, " cells [",percent,"%]\n", "(minFeatures=",minFeatures, "; maxFeatures=", maxFeatures, "; minCounts=" ,minCounts,  "; maxCounts=" ,maxCounts , "; maxMito=" ,maxMito, ") for ", unique(seuratObject$sample), "\n" )
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



#' Draw a FeatureHeatmap
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

## prepare data for cluster t.test from the deg list and do a cluster t-test
do_cluster_t_test <- function(seurat_subset, degs, group="condition", cluster="seurat_clusters"){
  gene_names<- names(table(rownames(degs)))
  #print(head(gene_names))
  p_values <- vector("list",length(gene_names))
  names(p_values) <- gene_names
  #gene_names <- row.names(cluster_subset)
  #if (celltype=="Adipocytes"){
  #  seurat_subset <- seurat_subset[,!seurat_subset$orig.ident=="D7"]
  #}
  group <- seurat_subset[[group]][,1]
  cluster <- seurat_subset[[cluster]][,1]
  for (gene in gene_names){
    y <- c(t(as.matrix(seurat_subset@assays$RNA[gene,])))
    test_info <- my.t.test.cluster(y, cluster, group)
    p_values[[gene]] <- test_info[nrow(test_info)]
  }
  return(p_values)
}

## added line 54-56 so that each group is tested if 
## only one observation is present and throw an error
my.t.test.cluster <- function (y, cluster, group, conf.int = 0.95) 
{
  group <- as.factor(group)
  cluster <- as.factor(cluster)
  s <- !(is.na(y) | is.na(cluster) | is.na(group))
  y <- y[s]
  cluster <- cluster[s]
  group <- group[s]
  n <- length(y)
  if (n < 2) 
    stop("n<2")
  gr <- levels(group)
  if (length(gr) != 2) 
    stop("must have exactly two treatment groups")
  n <- table(group)
  nc <- tapply(cluster, group, function(x) length(unique(x)))
  bar <- tapply(y, group, mean)
  u <- unclass(group)
  y1 <- y[u == 1]
  y2 <- y[u == 2]
  c1 <- factor(cluster[u == 1])
  c2 <- factor(cluster[u == 2])
  b1 <- tapply(y1, c1, mean)
  b2 <- tapply(y2, c2, mean)
  m1 <- table(c1)
  m2 <- table(c2)
  if (any(names(m1) != names(b1)))
    stop("logic error 1")
  if (any(names(m2) != names(b2)))
    stop("logic error 2")
  if (any(m2 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m2[m2 < 2]), collapse = " ")))
  if (any(m1 < 2))
    stop(paste("The following clusters contain only one observation:",
               paste(names(m1[m1 < 2]), collapse = " ")))
  M1 <- mean(y1)
  M2 <- mean(y2)
  ssc1 <- sum(m1 * ((b1 - M1)^2))
  ssc2 <- sum(m2 * ((b2 - M2)^2))
  if (nc[1] != length(m1))
    stop("logic error 3")
  if (nc[2] != length(m2))
    stop("logic error 4")
  df.msc <- sum(nc) - 2
  msc <- (ssc1 + ssc2)/df.msc
  v1 <- tapply(y1, c1, var)
  v2 <- tapply(y2, c2, var)
  ssw1 <- sum((m1 - 1) * v1)
  ssw2 <- sum((m2 - 1) * v2)
  df.mse <- sum(n) - sum(nc)
  mse <- (ssw1 + ssw2)/df.mse
  na <- (sum(n) - (sum(m1^2)/n[1] + sum(m2^2)/n[2]))/(sum(nc) - 
                                                        1)
  rho <- (msc - mse)/(msc + (na - 1) * mse)
  r <- max(rho, 0)
  C1 <- sum(m1 * (1 + (m1 - 1) * r))/n[1]
  C2 <- sum(m2 * (1 + (m2 - 1) * r))/n[2]
  v <- mse * (C1/n[1] + C2/n[2])
  v.unadj <- mse * (1/n[1] + 1/n[2])
  de <- v/v.unadj
  dif <- diff(bar)
  se <- sqrt(v)
  zcrit <- qnorm((1 + conf.int)/2)
  cl <- c(dif - zcrit * se, dif + zcrit * se)
  z <- dif/se
  P <- 2 * pnorm(-abs(z))
  stats <- matrix(NA, nrow = 20, ncol = 2, dimnames = list(c("N", 
                                                             "Clusters", "Mean", "SS among clusters within groups", 
                                                             "SS within clusters within groups", "MS among clusters within groups", 
                                                             "d.f.", "MS within clusters within groups", "d.f.", "Na", 
                                                             "Intracluster correlation", "Variance Correction Factor", 
                                                             "Variance of effect", "Variance without cluster adjustment", 
                                                             "Design Effect", "Effect (Difference in Means)", "S.E. of Effect", 
                                                             paste(format(conf.int), "Confidence limits"), "Z Statistic", 
                                                             "2-sided P Value"), gr))
  stats[1, ] <- n
  stats[2, ] <- nc
  stats[3, ] <- bar
  stats[4, ] <- c(ssc1, ssc2)
  stats[5, ] <- c(ssw1, ssw2)
  stats[6, 1] <- msc
  stats[7, 1] <- df.msc
  stats[8, 1] <- mse
  stats[9, 1] <- df.mse
  stats[10, 1] <- na
  stats[11, 1] <- rho
  stats[12, ] <- c(C1, C2)
  stats[13, 1] <- v
  stats[14, 1] <- v.unadj
  stats[15, 1] <- de
  stats[16, 1] <- dif
  stats[17, 1] <- se
  stats[18, ] <- cl
  stats[19, 1] <- z
  stats[20, 1] <- P
  attr(stats, "class") <- "t.test.cluster"
  stats
}

#Do Mean SEM Graph
#Autor Mariano

doMeanSemBarplot <- function(SeuratObject, Features, ListTest=NULL){ 
  require(ggplot2)
  require(ggpubr)
  require(tidyverse)
  #SEM function defintion
  SEM <- function(x) sqrt(var(x)/length(x))
  #create data frame with conditions from provided SeuratObject, aswell as original identifier of samples
  df<-data.frame(condition=SeuratObject$condition
                 ,orig.ident = SeuratObject$orig.ident)
  #get expression values for genes from individual cells, add to df
  for(i in Features){
    df[,i] <- SeuratObject@assays$RNA@data[i,]
  }
  
  #melt results 
  df.melt <- melt(df)
  #group results and summarize, also add/use SEM 
  df.melt.sum <- df.melt %>% 
    group_by(condition, variable) %>% 
    summarise(Mean = mean(value), SEM = SEM(value))
  #second dataframe containing mean values for individual samples
  df.melt.orig <- df.melt %>% 
    group_by(condition, variable, orig.ident) %>% 
    summarise(Mean = mean(value))
  
  #create comparison list for t.test, always against control, so please check your sample ordering
  # ,alternative add your own list as argument
  if (is.null(ListTest)) {
    #if ListTest is empty, so grep the ctrl conditions out of the list 
    # and define ListTest comparing every other condition with that ctrl condition
    cat("ListTest empty, comparing every sample with provided control")
    conditions <- unique(SeuratObject$condition)
    ctrl.condition <- conditions[grep(pattern = paste(c("CTRL","Ctrl"),collapse ="|")
                                      ,conditions)[1]]
    df.melt.sum$condition <- factor(df.melt.sum$condition
                                    ,levels = c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))
    #create ListTest
    ListTest <- list()
    for (i in 1:length(conditions)) {
      cndtn <- conditions[i] 
      if(cndtn!=ctrl.condition)
      {
        ListTest[[i]] <- c(ctrl.condition,cndtn)
      }
    }
  }
  #delete Null values, created by count index
  ListTest <- ListTest[!sapply(ListTest, is.null)]
  #create barplot with significance
  ggplot(df.melt.sum, aes(x = condition, y = Mean, fill = condition))+
    geom_col(color = "black")+
    geom_errorbar(aes(ymin = Mean-SEM, ymax = Mean+SEM), width = .2)+
    geom_jitter(data = df.melt.orig, aes(x=condition,y=Mean), size = 1, shape=1)+
    #ordering, control always first
    scale_x_discrete(limits=c(ctrl.condition,levels(factor(df.melt.sum$condition))[!(levels(factor(df.melt.sum$condition)) %in% ctrl.condition)]))+
    #t-test, always against control, using means from orig sample identifier
    stat_compare_means(data=df.melt.orig, comparisons = ListTest, method = "t.test", size=3)+
    facet_wrap(~variable, ncol = 9, scales = "free") +
    scale_fill_manual(values = c("#FFFFFF","#BDD7E7" ,"#6BAED6", "#3182BD", "#08519C","#0099CC","#00BFFF","#63D1F4","#00688B")
                      , name = "Condition")+
    labs( title = "", y = "Mean UMI") +
    theme_classic() +
    theme(axis.text.x = element_text(face = "bold", color = "black",angle = 45,hjust = 1, size = 10),
          axis.text.y = element_text(color = "black"),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 14, color = "black"),
          plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
          axis.line = element_line(color = "black", size = .6),
          legend.position = "bottom")
  
}

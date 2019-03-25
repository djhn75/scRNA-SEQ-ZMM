library(Seurat)

#' 
#' @author David John
#' @param seuratObject 
#' @return filtered seurat object
FilterDeadCells <- function(seuratObject){
  # The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
  # For non-UMI data, nCount_RNA represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
  # We use raw count data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
  mito.features <- grep(pattern = "^MT-", x = rownames(x = seuratObject), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = seuratObject, slot = 'counts'))
  seuratObject[['percent.mito']] <- percent.mito
  VlnPlot(object = seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
  
  seuratObject <- subset(x = seuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.05)
  return(seuratObject)
}

#' Import Single cell sequencing experiments into Seurat and perform normalisation and scale Data 
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE, performNormalisation=TRUE, performVariableGeneDetection=TRUE) {
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  
  seuratObject =CreateSeuratObject(counts = Matrix, min.cells= 3, min.features = 200)
  seuratObject<-RenameCells(object = seuratObject, add.cell.id = id)
  if (performNormalisation==TRUE) {
    seuratObject<-NormalizeData(object = seuratObject)
  }
  if(performVariableGeneDetection){
    seuratObject<-FindVariableFeatures(object = seuratObject, do.plot = FALSE, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  #seuratObject<-ScaleData(object = seuratObject)
  cat("Imported ", length(seuratObject@meta.data$orig.ident), " cells from ", pathway, "with ID ", id, "\n")
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



library(scImpute)
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
    scimpute(count_path = path.Matrix, out_dir = path.Imputed.Matrix, Kcluster = cluster, ncores=ncores, drop_thre = drop_thre)
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

library(Seurat)


#' Import Single cell sequencing experiments into Seurat and ronnarmalisation and scale Data 
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
Importer <- function(pathway,id, TenX=TRUE) {
  if (TenX) {
    Matrix <- Read10X(pathway)
  }  else{
    Matrix <- read.table(pathway,header = TRUE,sep = ",", dec = ".", row.names = 1)
  }
  
  seuratObject =CreateSeuratObject(raw.data = Matrix, min.cells= 3, min.genes = 200)
  seuratObject<-RenameCells(object = seuratObject, add.cell.id = id)
  seuratObject<-NormalizeData(object = seuratObject)
  seuratObject<-ScaleData(object = seuratObject)
  seuratObject<-FindVariableGenes(object = seuratObject, do.plot = FALSE)
  cat("Imported ", length(seuratObject@ident), " cells from ", pathway, "with ID ", id, "\n")
  return(seuratObject)
}

#' Import and combine several Single cell sequencing experiments into Seurat
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineSeuratObjectsCCA <- function(pathways,ids){
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  hvg<-c() #to save high variable genes
  seurat.objects<-list() #list to dave Seurat objects
  for (i in 1:length(pathways)) {
    seuratObject<-Importer(pathways[i],ids[i])
    hvg<-c(hvg,head(rownames(seuratObject@hvg.info), 1000)) #save variable genes
    seurat.objects<-c(seurat.objects,seuratObject)
  }
  hvg <- names(which(table(hvg) > 1))
  #check if the variable genes are expressed in all sample
  for (i in 1:length(seurat.objects)) {
    hvg <- hvg[hvg %in% rownames(seurat.objects[[i]]@scale.data)]
  }
  seuratObject<-RunMultiCCA(object.list = seurat.objects,genes.use = hvg,)
  cat("Merged Seurat object contains ", length(seuratObject@ident)," cells\n")
  return(seuratObject)
}


#' Import and combine several Single cell sequencing experiments into Seurat
#' @author David John
#' @param pathways A vector of pathways to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineSeuratObjects <- function(pathways,ids){
  Importer <- function(pathway,id) {
    Ten_X <- Read10X(pathway)
    seuratObject =CreateSeuratObject(raw.data = Ten_X, min.cells= 3, min.genes = 200)
    seuratObject<-RenameCells(object = seuratObject, add.cell.id = id)
    cat("Imported ", length(seuratObject@ident), " cells from ", pathway, "with ID ", id, "\n")
    return(seuratObject)
  }
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  for (i in 1:length(pathways)) {
    if (i<2) { 
      seuratObject1<-Importer(pathways[i],ids[i])
      next
    }
    seuratObject2<-Importer(pathways[i],ids[i])
    seuratObject1 <- MergeSeurat(object1 = seuratObject1, object2 = seuratObject2) 
  }
  cat("Merged Seurat object contains ", length(seuratObject1@ident)," cells\n")
  return(seuratObject1)
}



library(scImpute)
library(Seurat)
#' Read 10X results and impute the dataset with scImpute
#' @author David John
#' @param pathways Pathway to the cellrancer count output folder (contains barcodes.tsv, genes.tsv, matrix.mtx)
#' @param ids String that are assigned to the output matrix
imputeData <- function(pathway,id, cluster=12, ncores=20, drop_thre=0.5){
  for (i in 1:length(pathways)) {
    path.Matrix<-paste(pathways[i],"Matrix.csv",sep="/")
    path.Imputed.Matrix <- paste(pathways[i], "scImpute", ids[i], sep="/")
    Ten_X <- Read10X(pathways[i])
    write.csv(as.data.frame(as.matrix(Ten_X)), file = path.Matrix)
    scimpute(count_path = path.Matrix,
           out_dir = path.Imputed.Matrix, Kcluster = cluster, ncores=ncores, drop_thre = drop_thre)
    cat("Wrote imputed Martix to ", path.Imputed.Matrix)
  }
}



#' Import and combine imputet single cell Matrices into Seurat
#' @author David John
#' @param pathways A vector of pathways to scimputed count marices 
#' @param ids Vector of strings that are assigned to the concordant cells
#' @return Merged seurat object
combineScImputedSeuratObjectsCCA <- function(pathways,ids){
  
  if (length(pathways)!=length(ids)) {stop(" pathways and ids vector need to have the same length")  }
  hvg<-c() #to save high variable genes
  seurat.objects<-list() #list to dave Seurat objects
  #loop through directories of pathways
  for (i in 1:length(pathways)) {
    matrixPath<-paste(pathways[i],"scImpute",paste(ids[i],"scimpute_count.csv",sep=""),sep = "/")
    seuratObject<-Importer(matrixPath,ids[i], TenX = FALSE)
    hvg<-c(hvg,head(rownames(seuratObject@hvg.info), 1000)) #save variable genes
    seurat.objects<-c(seurat.objects,seuratObject)
  }
  hvg <- names(which(table(hvg) > 1))
  #check if the variable genes are expressed in all sample
  for (i in 1:length(seurat.objects)) {
    hvg <- hvg[hvg %in% rownames(seurat.objects[[i]]@scale.data)]
  }
  
  seuratObject<-RunMultiCCA(object.list = seurat.objects,genes.use = hvg,)
  cat("Merged Seurat object contains ", length(seuratObject@ident)," cells\n")
  return(seuratObject)
}

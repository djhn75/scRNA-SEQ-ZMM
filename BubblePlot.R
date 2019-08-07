Data = AMI.Simone
library(hash)
Genes <- c("Pdgfra", "Col4a1", "Adamts5", "Lamb1", "Dpep1", "Ms4a4d", "Medag", "Col3a1", "Mmp2", "Dpt", "Dcn", "Dkk3", "Comp", "Tbx20",
           "Meox1", "Prg4", "Frzb", "Tnc", "Col1a1", "Lox", "Postn", "Acta2","Lyz2", "Lgals3", "Mrc1", "Cd68", "Cd14", "Msr1", "Cd209a",
           "Cd83", "Napsa", "Cd74", "Ly6g", "Cxcr2", "Prox1", "Tie1", "Lyve1", "Cldn5", "Vwf", "Eng", "Emcn", "Fabp4", "Pecam1", "Cdh5",
           "Ccl5", "Nkg7", "Ptprc", "Klrc1", "Ctla4", "Klre1", "Cd3g", "Trdc", "Icos", "Cd3e", "Lat", "Lef1", "Tcf7", "Iglc1", "Pax5",
           "Iglc2", "Iglc3", "Cd79a", "Cd79b", "Cd19", "Plp1", "Kcna1", "Kcna2", "Cd59a", "Rgs5", "Tagln", "Myh11", "Vtn", "Notch3",
           "Pdgfrb", "Cspg4", "Des", "Pln", "Nkain4", "Krt8", "Krt19", "Krt18")

Annotations <- c(rep(0,22), rep(1,12), rep(2,10), rep(3,6), rep(4,7), rep(5,7), rep(6,4), rep(7,9), rep(8,4))
Celltypes <- c("Fibroblasts", "Fibroblasts", "NK-cells", "Monocytes", "Endothelial cells", "Monocytes", "T-cells", "Fibroblasts", "B-cells", "Fibroblasts", "Monocytes", "Endothelial cells", "T-cells", "Schwann cells", "Pericytes/Smooth muscle cells", "Pericardial cells", "Monocytes")

for (i in Genes){
  if(length(Data@data[i,] == length(rownames(Data@meta.data)))){print(i)}
}


Markers.TSNE <- FindAllMarkers(Data, genes.use = Genes, logfc.threshold = 0, test.use = "bimod", min.pct = 0)

A <- data.frame(Gene = Markers.TSNE$gene, log_FC = Markers.TSNE$avg_logFC, 
                p.value = Markers.TSNE$p_val_adj, Cluster = Markers.TSNE$cluster)
# Replace 0 with 1.0e-300
for (i in rownames(A)){
  if(A[i,"p.value"] == 0){
    A[i,"p.value"] <- 1.0e-300
  }
}

gene.cluster <- hash(keys = Genes, values = Annotations)
cells.cluster <- hash(keys = as.character(c(0:16)), values = Celltypes)

y = NULL
x = NULL
for(i in rownames(A)){
  y <- c(y, gene.cluster[[as.character(A[i,"Gene"])]])
  x <- c(x, cells.cluster[[as.character(A[i,"Cluster"])]])
}

A <- data.frame(A, log.10.p = -log10(A$p.value), Cluster.Gene = y, Cluster.Celltype = x)

sort<- A[order(-A$log_FC), ]

A$Gene <- factor(A$Gene, levels = unique(sort$Gene))

ggplot(A, aes(x=Cluster, y= Gene, color = log_FC, size = log.10.p))+ 
  geom_point(pch = 19)+ 
  scale_color_gradient2(low = "blue1", mid = "grey80", high = "red1", name = "ln(Fold change)")+
  scale_size(name = "log10(p-value)")+ 
  ggtitle("Marker in all cells")+
  facet_grid(Cluster.Celltype~Cluster.Gene, space = "free", scales = "free")+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, face = "italic"), legend.position = "top")
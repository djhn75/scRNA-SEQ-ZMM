loadRDS("/media/ATLAS_NGS_storage/Wesley/CircRes2019/CircResBothCohortsMonocytes.final.Rds")
library(Seurat)
FindMarkers()

packageVersion("Seurat")
CircResBothCohortsMonocytes.final<-readRDS("/media/ATLAS_NGS_storage/Wesley/CircRes2019/CircResBothCohortsMonocytes.final.Rds")




StudyTable<-read.table("/media/ATLAS_NGS_storage/Wesley/Wesley-Study-Table-for statistics_08052019_updated_09.09.19.csv", skip = 1, header = TRUE, sep = "\t")
rownames(StudyTable)<-StudyTable$ScRNA.seq.ID

CircResBothCohortsMonocytes.final<-SetAllIdent(CircResBothCohortsMonocytes.final, id = "newSamplenames")

meta.data<-c()
colname<-"Diabetes.Mellitus"

addStudyTable <- function(seuratObject, studyTable, colname){
  meta.data<-c()
    for (ident in seuratObject@ident) {
      meta.data<-c(meta.data,as.character(studyTable[ident, colname]))
    
  }
  names(meta.data)<-seuratObject@cell.names
  seuratObject<-AddMetaData(object = seuratObject,metadata = meta.data,col.name = colname);rm(meta.data)
  
  return(seuratObject)
}

variablesToAdd<-c("Intervention","leading.diagnosis","Diabetes.Mellitus","Insulin",
                  "Renal.Insufficiency","Hypertension","Hyperlipidemia","Smoking",
                  "s.p.myocardial.infarction","s.p.PCI","s.p.CABG.Bypass.surgery",
                  "s.p.apoplex.stroke","s.p.pAVK","s.p.apoplex.stroke","s.p.pAVK",
                  "s.p.PTA","AGE","Gender","ICM","Ethiology..ICM.vs.non.ICM.","LVEF",
                  "LVEF..","NYHA","Systolic.BP","Weight..kg.","LBBB.Left.Bundle.Branch.Block..QRS.120ms.","Atrial.fibrillation")

for (tmp in colnames(StudyTable)){
  CircResBothCohortsMonocytes.final<-addStudyTable(CircResBothCohortsMonocytes.final, studyTable = StudyTable, colname = tmp)
}

saveRDS(CircResBothCohortsMonocytes.final,"/media/ATLAS_NGS_storage/Wesley/CircRes2019/CircResBothCohortsMonocytes.final.metaData.Rds")
CircResBothCohortsMonocytes.final<-readRDS("/media/ATLAS_NGS_storage/Wesley/CircRes2019/CircResBothCohortsMonocytes.final.metaData.Rds")

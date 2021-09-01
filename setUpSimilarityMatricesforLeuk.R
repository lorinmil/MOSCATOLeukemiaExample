
#---------------Specify inputs--------------------#
library(Seurat)
library(patchwork)
library(rTensor)

#Source the function to process the files
source("function_seuratToMatrix.R")

outputLoc <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/"





#---------------Read in necessary files-----------#
immune.combined <- readRDS(paste0(outputLoc, "filesFromMerge/seurat_combined_all.rda"))
metaCombined <- immune.combined@meta.data
rm(immune.combined)

#Clean out the rownames
cellNames <- rownames(metaCombined)
healthyIdx <- which(substr(cellNames, 1, 8)=="Healthy_")
leukIdx <- which(substr(cellNames, 1, 8)!="Healthy_")
cellNames[healthyIdx] <- gsub("^(Healthy__)?", "", cellNames[healthyIdx])
cellNames[leukIdx] <- gsub("^(leuk__)?", "", cellNames[leukIdx])
rownames(metaCombined) <- cellNames




#-------Convert from Seurat to R data objects-----#
seuratInfo <- seuratToMatrix(seuratLoc_="C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/SeuratObjects",
                             outputLoc_="C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/rObjects",
                             integratedSeuratMeta_=metaCombined)
saveRDS(seuratInfo, "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/rObjects/seuratInfo.rda")




#-------Find which chromosome the RNA live-----#
library('biomaRt')
mart <- useEnsembl(biomart = "ensembl",
                   dataset = "hsapiens_gene_ensembl")
geneInfo <- getBM(attributes = c("external_gene_name", "chromosome_name", "description"),
                  filters = "external_gene_name",
                  values = seuratInfo$commonRNA,
                  mart = mart)
geneInfo <- geneInfo[geneInfo$chromosome_name %in% 1:22,]
table(geneInfo$chromosome_name[geneInfo$chromosome_name %in% 1:22])




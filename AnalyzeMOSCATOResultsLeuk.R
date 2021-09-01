


folderLoc_ <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication"
outLoc_ <- paste0(folderLoc_, "/MOSCATO/output/v2")
geneInfo <- readRDS("C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/rObjects/seuratInfo.rda")

#This is the cluster order:
# [1] "similarity_clust_0.rda"  "similarity_clust_1.rda"  "similarity_clust_14.rda" "similarity_clust_2.rda"  "similarity_clust_3.rda" 
# [6] "similarity_clust_4.rda"  "similarity_clust_5.rda"  "similarity_clust_6.rda"




#------------Cluster 0--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_1.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)




#------------Cluster 1--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_2.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)





#------------Cluster 14--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_3.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)






#------------Cluster 2--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_4.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)





#------------Cluster 3--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_5.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)






#-----Cluster 4 --------------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_6.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)






#------------Cluster 5--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_7.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)






#------------Cluster 6--------------------#
networkResults <- readRDS(paste0(outLoc_, "/scNetworkResults_8.rda"))
selectRNA <- geneInfo$commonRNA[networkResults$finalBeta[[1]]!=0]
selectADT <- geneInfo$commonADT[networkResults$finalBeta[[2]]!=0]
View(selectRNA)
View(selectADT)

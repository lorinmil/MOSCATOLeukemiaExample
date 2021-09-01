
library(Seurat)
library(rTensor)

seuratToMatrix <- function(seuratLoc_, 
                           seuratFiles_=NULL,
                           outputLoc_,
                           integratedSeurat_=NULL,
                           integratedSeuratMeta_=NULL)
{
  
  #List the files from the seurat location
  if(!is.null(seuratFiles_)){
    seuratFiles <- paste0(seuratLoc_, "/", seuratFiles)
  }else {
    seuratFiles <- paste0(seuratLoc_, "/", list.files(seuratLoc_))
  }
  
  #Get the Cell type info
  if(!is.null(integratedSeurat_)){
    cellInfo <- integratedSeurat_@meta.data
  }else if(!is.null(integratedSeuratMeta_)){
    cellInfo <- integratedSeuratMeta_
  }else {
    cat("Need to input either integratedSeurat_ or integratedSeuratMeta_")
  }
  
  
  
  #STEP1: Loop through each of the Seurat files to see what they have in common
  for(file in 1:length(seuratFiles)){
    
    #Read the Seurat file
    seurat <- readRDS(seuratFiles[file])
    
    #Which cells does this individual belong to?
    cellsFromFile <- which(rownames(cellInfo) %in% colnames(seurat@assays[[1]]))
    #What are the cell cluster assignments?
    cellCluster <- cellInfo[cellsFromFile, "seurat_clusters"]
    
    
    #Add to the list
    if(file==1){
      commonCellClusts <- unique(cellCluster)
      commonRNA <- rownames(seurat@assays[["RNA"]])
      commonADT <- rownames(seurat@assays[["ADT"]])
    }else {
      commonCellClusts <- unique(intersect(commonCellClusts, cellCluster))
      commonRNA <- unique(intersect(commonRNA, rownames(seurat@assays[["RNA"]])))
      commonADT <- unique(intersect(commonADT, rownames(seurat@assays[["ADT"]])))
    }
    
    cat("STEP1: Finished with Seurat Object ", file, "\n")
  } #Loop to the next Seurat file
  
  
  
  
  #STEP2: Loop through each of the cell clusters and create a similarity tensor
  for(clust in commonCellClusts){
    
    #Initialize
    similarityList <- list()
    
    for(file in 1:length(seuratFiles)){
      
      #Read the Seurat file
      seurat <- readRDS(seuratFiles[file])
      
      #Which cells does this individual belong to?
      cellsFromFile <- which(rownames(cellInfo) %in% colnames(seurat@assays[[1]]))
      #What are the cell cluster assignments?
      cellCluster <- cellInfo[cellsFromFile, "seurat_clusters"]
      
      #Pull out the cells for this cluster
      X_k <- t(as.matrix(seurat@assays[["RNA"]][commonRNA, as.character(cellCluster) %in% clust]))
      G_k <- t(as.matrix(seurat@assays[["ADT"]][commonADT, as.character(cellCluster) %in% clust]))
      
      #Calculate the similarity
      similarityList[[file]] <- cor(x=X_k, y=G_k, method="pearson")#, use="pairwise.complete.obs")
      #Convert any NAs to 0s
      similarityList[[file]][is.na(similarityList[[file]])] <- 0
    } #Loop to the next seurat file
    
    p <- ncol(X_k)
    q <- ncol(G_k)
    
    #Put the similarity matrix list into a tensor
    similarityArray <- array(data=unlist(similarityList), dim=c(p,q,length(seuratFiles)))
    X_tensor <- rTensor::as.tensor(similarityArray)
    
    #Export this cell cluster
    saveRDS(X_tensor, paste0(outputLoc_, "/similarity_clust_", clust, ".rda"))
    
    cat("STEP2: Finished with cluster ", clust, "\n")
    
  } #Loop to the next cluster
  
  #Return which clusters and features were used
  return(list(commonCellClusts=commonCellClusts,
              commonRNA=commonRNA,
              commonADT=commonADT))
  
  
} #end seuratToMatrix
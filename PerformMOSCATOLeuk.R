



folderLoc_ <- "/user/lorinmil/realDataApplication"

#Set up files and folders to load
outLoc_ <- paste0(folderLoc_, "/v2/output")
dataLoc_ <- paste0(folderLoc_, "/data")
similarityFiles <- list.files(dataLoc_)
similarityFiles <- paste0(dataLoc_, "/", similarityFiles)


###Import libraries###
library(ggplot2)
library(rTensor, lib.loc="/user/lorinmil/rlibs")
library(glmnet, lib.loc="/user/lorinmil/rlibs")
library(caret, lib.loc="/user/lorinmil/rlibs")
library(parallel)


###SOURCE SUPPORTING FUNCTIONS###
source("supportingMOSCATOFunctions.R")
source("tuneMOSCATO.R")


for(file in 1:length(similarityFiles)){
  
  time1 <- Sys.time()
  tuned <- tuneMOSCATO(Z_tensor=readRDS(similarityFiles[file]),
                       y=c(rep(0, 14), rep(1, 7)),
                       distY="binomial",
                       #Tensor regression parameters
                       epsilon=0.1,
                       minIter=5,
                       maxIter=25,
                       phi=0.001,
                       alphaGrid=c(0.01,0.05,0.1),
                       minSelect=c(10, 1),
                       numIterForStability=50, #Number of iterations to establish stability
                       sizeOfSubsampleForStability=20 #How many subjects should be used for the subsampling in StARS?
  )
  testNetwork <- tensorElasticNet(Z_tensor=readRDS(similarityFiles[file]),
                                  y=c(rep(0, 14), rep(1, 7)),
                                  maxIter=20,
                                  glmnetAlpha=tuned$tunedAlpha,
                                  maxSelect=tuned$tunedMax
  )
  time2 <- Sys.time()
  #Add the runtime
  timeDiff <- time2 - time1
  testNetwork$runtime <- timeDiff
  
  #Save the network results
  saveResults <- list(tuned=tuned,
                      MOSCATO=testNetwork)
  saveRDS(saveResults, paste0(outLoc_, "/scNetworkResults_", file, ".rda"))
  
}
  

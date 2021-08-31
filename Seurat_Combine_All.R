

#Match cell types across projects
library(Seurat)
library(patchwork)

outputLoc <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/"

# memory.limit(size=54000)


#--------------------------------------------------------------#
#-----------------STEP1: Setup the Seurat objects--------------#
#--------------------------------------------------------------#

#What are all the Seurat files?
seuratObjs <- list.files(paste0(outputLoc, "SeuratObjects"))

ifnb.list <- list()
for(i in 1:length(seuratObjs)){
  ifnb.list[[i]] <- readRDS(paste0(outputLoc, "SeuratObjects/", seuratObjs[i]))
}


# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})

#Rename the cells to have healthy vs leukemia
leukemia <- c(21, 15, 16, 17, 18, 19, 20)
healthy <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
for(idx in leukemia){
  ifnb.list[[idx]] <- RenameCells(object=ifnb.list[[idx]],
                                  add.cell.id="leuk_")
}
for(idx in healthy){
  ifnb.list[[idx]] <- RenameCells(object=ifnb.list[[idx]],
                                  add.cell.id="Healthy_")
}

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list,
                                      nfeatures=500)


#Used since now using PCA instead of CCA
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list,
                                      nfeatures=500)

#Save
saveRDS(ifnb.list, paste0(outputLoc, "filesFromMerge/seurat_list_all.rda"))
saveRDS(features, paste0(outputLoc, "filesFromMerge/seurat_features_all.rda"))

#Read the objects (if the above has already been ran and returning to R)
ifnb.list <- readRDS(paste0(outputLoc, "filesFromMerge/seurat_list_all.rda"))
features <- readRDS(paste0(outputLoc, "filesFromMerge/seurat_features_all.rda"))

#Find the anchors for integration
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")

#Remove to save memory
rm(ifnb.list)
rm(features)
rm(healthy)
rm(leukemia)
rm(seuratObjs)
rm(idx)

#Read if returning to R
immune.anchors <- readRDS(paste0(outputLoc, "filesFromMerge/seurat_anchors_all.rda"))

#Create an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

#Save the objects
saveRDS(immune.anchors, paste0(outputLoc, "filesFromMerge/seurat_anchors_all.rda"))
saveRDS(immune.combined, paste0(outputLoc, "filesFromMerge/seurat_combined_all.rda"))




# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DimPlot(immune.combined, reduction = "umap", label = TRUE) 
# Add to metadata
immune.combined@meta.data$Outcome <- "Leukemia"
immune.combined@meta.data$Outcome[grepl("^(Healthy_)", rownames(immune.combined@meta.data))] <- "Healthy"
DimPlot(immune.combined, reduction = "umap", group.by = "Outcome")
DimPlot(immune.combined, reduction = "umap", split.by = "Outcome", label = TRUE)

#Save after re-processing the integration
saveRDS(immune.combined, paste0(outputLoc, "filesFromMerge/seurat_combined_all.rda"))

#Read if returning to R
immune.combined <- readRDS(paste0(outputLoc, "filesFromMerge/seurat_combined_all.rda"))

#Find the important genes for the interesting cell types with only healthy subjects
adt_markers_healthy <- FindMarkers(immune.combined, ident.1 = 15, assay = "ADT")
rna_markers_healthy <- FindMarkers(immune.combined, ident.1 = 15, assay = "RNA")
rna_markers_healthy[1:15,]
#Make the plots for those genes
# p1 <- FeaturePlot(immune.combined_healthy, "adt_CD27", cols = c("lightgrey", "darkgreen")) + ggtitle("CD27 protein")
FeaturePlot(immune.combined, "rna_SLC8A1-AS1") + ggtitle("SLC8A1-AS1 RNA")
FeaturePlot(immune.combined, "rna_CMTM8") + ggtitle("CMTM8 RNA")
FeaturePlot(immune.combined, "rna_CMTM7") + ggtitle("CMTM7 RNA")
FeaturePlot(immune.combined, "rna_FHIT") + ggtitle("FHIT RNA")
FeaturePlot(immune.combined, "rna_MYLK") + ggtitle("MYLK RNA")
write.csv(rownames(rna_markers_healthy), paste0(outputLoc, "filesFromMerge/seurat_combined_all_findMarkers15.csv"))



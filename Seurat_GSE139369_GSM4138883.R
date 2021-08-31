
library(Seurat)
library(ggplot2)
library(patchwork)

#Specify where the data is stored
outputLoc <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/"

#------------------------------------------------------------------------#
#----------------STEP 1: Load Data---------------------------------------#
#------------------------------------------------------------------------#

#Import the data
study4_RNA <- readRDS(paste0(outputLoc, "GSE139369/GSM4138883_scRNA_MPAL4_T1.rds"))
study4_protein <- readRDS(paste0(outputLoc, "GSE139369/GSM4138883_scADT_MPAL4_T1.rds"))

#Transpose so columns are cells (for Seurat) and make a sparse object
study4_RNA <- as.sparse(study4_RNA)
study4_protein <- as.sparse(study4_protein)

#Make sure that the cells match across data types
all.equal(colnames(study4_RNA), colnames(study4_protein)) #False!
cellTags_RNA <- colnames(study4_RNA)
cellTags_protein <- colnames(study4_protein)
matchingCells <- cellTags_RNA[cellTags_RNA %in% cellTags_protein]
study4_RNA <- study4_RNA[, match(matchingCells, cellTags_RNA)]
study4_protein <- study4_protein[, match(matchingCells, cellTags_protein)]
all.equal(colnames(study4_RNA), colnames(study4_protein)) #Now it's good!




#------------------------------------------------------------------------#
#----------------STEP 2: Set up Seurat Object----------------------------#
#------------------------------------------------------------------------#

#Add RNA
study4 <- CreateSeuratObject(counts = study4_RNA)
Assays(study4) #Good! Says RNA which is the default

#Add proteins
adt_assay <- CreateAssayObject(counts = study4_protein)
study4[["ADT"]] <- adt_assay
Assays(study4) #Good! Now says RNA and ADT (which is the proteins)

#Extract a list of features measured in the ADT assay
rownames(study4[["ADT"]])
#List the current default assay
DefaultAssay(study4)




#------------------------------------------------------------------------#
#----------------STEP 3: Cluster Cells by RNA----------------------------#
#------------------------------------------------------------------------#

#NOTE: Only visualizing and normalizing the RNA data
study4 <- NormalizeData(study4)
study4 <- FindVariableFeatures(study4)
study4 <- ScaleData(study4)
study4 <- RunPCA(study4, verbose = FALSE)
study4 <- FindNeighbors(study4, dims = 1:30)
study4 <- FindClusters(study4, resolution = 0.8, verbose = FALSE)
study4 <- RunUMAP(study4, dims = 1:30)
DimPlot(study4, label = TRUE) #cluster 5 looks much different than the others. Maybe only 2 cell types?




#------------------------------------------------------------------------#
#----------------STEP 4: Visualize Multi-Omic Side-by-Side---------------#
#------------------------------------------------------------------------#

# Normalize ADT data,
DefaultAssay(study4) <- "ADT"
study4 <- NormalizeData(study4, normalization.method = "CLR", margin = 2)
DefaultAssay(study4) <- "RNA"

# Note that the following command is an alternative but returns the same result
# study4 <- NormalizeData(study4, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD19 levels for RNA and protein By setting the default assay, we can
# visualize one or the other (JUST AN EXAMPLE GENE!). Only CD14, CD19, and CD4 are in RNA
DefaultAssay(study4) <- "ADT"
p1 <- FeaturePlot(study4, "CD14", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(study4) <- "RNA"
p2 <- FeaturePlot(study4, "CD14") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2

# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(study4[["RNA"]])
Key(study4[["ADT"]])
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(study4, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(study4, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2




#------------------------------------------------------------------------#
#----------------STEP 5: Identify Proteins from RNA Clusters-------------#
#------------------------------------------------------------------------#

# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the
# surface
VlnPlot(study4, "adt_CD19") #Cluster 5 is showing up! We could see this visually from previous plots

# we can also identify alternative protein and RNA markers for this cluster through differential
# expression
adt_markers <- FindMarkers(study4, ident.1 = 5, assay = "ADT")
rna_markers <- FindMarkers(study4, ident.1 = 5, assay = "RNA")

head(adt_markers)
head(rna_markers)




#------------------------------------------------------------------------#
#----------------STEP 6: Additional Visualizations-----------------------#
#------------------------------------------------------------------------#

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(study4, feature1 = "adt_CD19", feature2 = "adt_CD3")

# view relationship between protein and RNA
FeatureScatter(study4, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(study4, feature1 = "adt_CD4", feature2 = "adt_CD8")




#------------------------------------------------------------------------#
#----------------STEP 7: Save the Seurat Object--------------------------#
#------------------------------------------------------------------------#

saveRDS(study4, paste0(outputLoc, "SeuratObjects/seuratObj_GSE139369_GSM4138883.rda"))


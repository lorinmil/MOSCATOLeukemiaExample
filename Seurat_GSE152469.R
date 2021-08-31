
library(Seurat)
library(ggplot2)
library(patchwork)

#Specify where the data is stored
outputLoc <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/"

#------------------------------------------------------------------------#
#----------------STEP 1: Load Data---------------------------------------#
#------------------------------------------------------------------------#

#Import the data
study3_RNA <- read.table(paste0(outputLoc, "GSE152469/GSM4616298_Exp_UMI_CLL_CPT1_M0.tsv"), sep = '\t', header = TRUE)
study3_protein <- read.csv(paste0(outputLoc, "GSE152469/GSM4616298_Ab_CLL_CPT1_M0.csv"), sep="")

#Pull out the cell info and make it a sparse object for Seurat
cellTags_RNA <- study3_RNA[,1]
cellTags_protein <- study3_protein[,1]

#Transpose so columns are cells (for Seurat) and make a sparse object
study3_RNA <- as.sparse(t(study3_RNA[,2:ncol(study3_RNA)]))
study3_protein <- as.sparse(t(study3_protein[,2:ncol(study3_protein)]))

#Make sure that the cells match across data types
cellTags_RNA_mod <- substr(cellTags_RNA, 1, nchar(cellTags_RNA)-2)
matchingCells <- cellTags_RNA_mod[cellTags_RNA_mod %in% cellTags_protein]
study3_RNA <- study3_RNA[, match(matchingCells, cellTags_RNA_mod)]
study3_protein <- study3_protein[, match(matchingCells, cellTags_protein)]

#Make the cell labels the column names
colnames(study3_RNA) <- matchingCells
colnames(study3_protein) <- matchingCells




#------------------------------------------------------------------------#
#----------------STEP 2: Set up Seurat Object----------------------------#
#------------------------------------------------------------------------#

#Add RNA
study3 <- CreateSeuratObject(counts = study3_RNA)
Assays(study3) #Good! Says RNA which is the default

#Add proteins
adt_assay <- CreateAssayObject(counts = study3_protein)
study3[["ADT"]] <- adt_assay
Assays(study3) #Good! Now says RNA and ADT (which is the proteins)

#Extract a list of features measured in the ADT assay
rownames(study3[["ADT"]])
#List the current default assay
DefaultAssay(study3)




#------------------------------------------------------------------------#
#----------------STEP 3: Cluster Cells by RNA----------------------------#
#------------------------------------------------------------------------#

#NOTE: Only visualizing and normalizing the RNA data
study3 <- NormalizeData(study3)
study3 <- FindVariableFeatures(study3)
study3 <- ScaleData(study3)
study3 <- RunPCA(study3, verbose = FALSE)
study3 <- FindNeighbors(study3, dims = 1:30)
study3 <- FindClusters(study3, resolution = 0.8, verbose = FALSE)
study3 <- RunUMAP(study3, dims = 1:30)
DimPlot(study3, label = TRUE) #cluster 5 looks much different than the others. Maybe only 2 cell types?




#------------------------------------------------------------------------#
#----------------STEP 4: Visualize Multi-Omic Side-by-Side---------------#
#------------------------------------------------------------------------#

# Normalize ADT data,
DefaultAssay(study3) <- "ADT"
study3 <- NormalizeData(study3, normalization.method = "CLR", margin = 2)
DefaultAssay(study3) <- "RNA"

# Note that the following command is an alternative but returns the same result
study3 <- NormalizeData(study3, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD19 levels for RNA and protein By setting the default assay, we can
# visualize one or the other (JUST AN EXAMPLE GENE!). Only CD14, CD19, and CD4 are in RNA
DefaultAssay(study3) <- "ADT"
p1 <- FeaturePlot(study3, "CD56", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(study3) <- "RNA"
p2 <- FeaturePlot(study3, "CD56") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2

# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(study3[["RNA"]])
Key(study3[["ADT"]])
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(study3, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(study3, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2




#------------------------------------------------------------------------#
#----------------STEP 5: Identify Proteins from RNA Clusters-------------#
#------------------------------------------------------------------------#

# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the
# surface
VlnPlot(study3, "adt_CD19") #Cluster 5 is showing up! We could see this visually from previous plots

# we can also identify alternative protein and RNA markers for this cluster through differential
# expression
adt_markers <- FindMarkers(study3, ident.1 = 5, assay = "ADT")
rna_markers <- FindMarkers(study3, ident.1 = 5, assay = "RNA")

head(adt_markers)
head(rna_markers)




#------------------------------------------------------------------------#
#----------------STEP 6: Additional Visualizations-----------------------#
#------------------------------------------------------------------------#

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(study3, feature1 = "adt_CD19", feature2 = "adt_CD3")

# view relationship between protein and RNA
FeatureScatter(study3, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(study3, feature1 = "adt_CD4", feature2 = "adt_CD8")

# Let's look at the raw (non-normalized) ADT counts. You can see the values are quite high,
# particularly in comparison to RNA values. This is due to the significantly higher protein copy
# number in cells, which significantly reduces 'drop-out' in ADT data
FeatureScatter(study3, feature1 = "adt_CD4", feature2 = "adt_CD8", slot = "counts")




#------------------------------------------------------------------------#
#----------------STEP 7: Save the Seurat Object--------------------------#
#------------------------------------------------------------------------#

saveRDS(study3, paste0(outputLoc, "SeuratObjects/seuratObj_GSE152469.rda"))


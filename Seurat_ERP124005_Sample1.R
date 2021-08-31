
library(Seurat)
library(ggplot2)
library(patchwork)
library('biomaRt')

#Specify where the data is stored
outputLoc <- "C:/Users/lorin/Box/Thesis_Lorin/realDataApplication/data/"

#------------------------------------------------------------------------#
#----------------STEP 1: Load Data---------------------------------------#
#------------------------------------------------------------------------#

#Read the cell annotations
study2_cellAnnotations <- read.csv(paste0(outputLoc, "ERP124005/CZI.PBMC.cell.annotations.csv"))
donors <- unique(study2_cellAnnotations$Donor_of_Origin)

#Import the data. Only want baseline data
study2_RNA <- readRDS(paste0(outputLoc, "ERP124005/CZI.PBMC.RNA.matrix.Rds"))
study2_RNA <- study2_RNA[, (colnames(study2_RNA) %in% study2_cellAnnotations$X[study2_cellAnnotations$HTO_Barcodes=="Baseline"]) & 
                           (colnames(study2_RNA) %in% study2_cellAnnotations$X[study2_cellAnnotations$Donor_of_Origin==donors[1]])]

study2_protein <- readRDS(paste0(outputLoc, "ERP124005/CZI.PBMC.ADT.matrix.Rds"))
study2_protein <- study2_protein[, (colnames(study2_protein) %in% study2_cellAnnotations$X[study2_cellAnnotations$HTO_Barcodes=="Baseline"]) & 
                                   (colnames(study2_protein) %in% study2_cellAnnotations$X[study2_cellAnnotations$Donor_of_Origin==donors[1]])]

#Transpose so columns are cells (for Seurat) and make a sparse object
study2_RNA <- as.sparse(study2_RNA)
study2_protein <- as.sparse(study2_protein)

#Clean up gene names
rownames(study2_protein) <- gsub("\\-(.*)", "", rownames(study2_protein))
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
geneInfo <- getBM(filters = "ensembl_gene_id", 
                  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                  values = rownames(study2_RNA), mart = mart)
study2_RNA <- study2_RNA[rownames(study2_RNA) %in% geneInfo$ensembl_gene_id,]
rownames(study2_RNA) <- geneInfo$external_gene_name[match(rownames(study2_RNA), geneInfo$ensembl_gene_id)]
study2_RNA <- study2_RNA[rownames(study2_RNA)!="",]

#Make sure that the cells match across data types
all.equal(colnames(study2_RNA), colnames(study2_protein)) #TRUE




#------------------------------------------------------------------------#
#----------------STEP 2: Set up Seurat Object----------------------------#
#------------------------------------------------------------------------#

#Add RNA
study2 <- CreateSeuratObject(counts = study2_RNA)
Assays(study2) #Good! Says RNA which is the default

#Add proteins
adt_assay <- CreateAssayObject(counts = study2_protein)
study2[["ADT"]] <- adt_assay
Assays(study2) #Good! Now says RNA and ADT (which is the proteins)

#Extract a list of features measured in the ADT assay
rownames(study2[["ADT"]])
#List the current default assay
DefaultAssay(study2)




#------------------------------------------------------------------------#
#----------------STEP 3: Cluster Cells by RNA----------------------------#
#------------------------------------------------------------------------#

#NOTE: Only visualizing and normalizing the RNA data
study2 <- NormalizeData(study2)
study2 <- FindVariableFeatures(study2)
study2 <- ScaleData(study2)
study2 <- RunPCA(study2, verbose = FALSE)
study2 <- FindNeighbors(study2, dims = 1:30)
study2 <- FindClusters(study2, resolution = 0.8, verbose = FALSE)
study2 <- RunUMAP(study2, dims = 1:30)
DimPlot(study2, label = TRUE) #cluster 5 looks much different than the others. Maybe only 2 cell types?




#------------------------------------------------------------------------#
#----------------STEP 4: Visualize Multi-Omic Side-by-Side---------------#
#------------------------------------------------------------------------#

# Normalize ADT data,
DefaultAssay(study2) <- "ADT"
study2 <- NormalizeData(study2, normalization.method = "CLR", margin = 2)
DefaultAssay(study2) <- "RNA"

# Note that the following command is an alternative but returns the same result
# study2 <- NormalizeData(study2, normalization.method = "CLR", margin = 2, assay = "ADT")

# Now, we will visualize CD19 levels for RNA and protein By setting the default assay, we can
# visualize one or the other (JUST AN EXAMPLE GENE!). Only CD14, CD19, and CD4 are in RNA
DefaultAssay(study2) <- "ADT"
p1 <- FeaturePlot(study2, "CD14", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(study2) <- "RNA"
p2 <- FeaturePlot(study2, "CD14") + ggtitle("CD19 RNA")

# place plots side-by-side
p1 | p2

# Alternately, we can use specific assay keys to specify a specific modality Identify the key
# for the RNA and protein assays
Key(study2[["RNA"]])
Key(study2[["ADT"]])
# Now, we can include the key in the feature name, which overrides the default assay
p1 <- FeaturePlot(study2, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
p2 <- FeaturePlot(study2, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2




#------------------------------------------------------------------------#
#----------------STEP 5: Identify Proteins from RNA Clusters-------------#
#------------------------------------------------------------------------#

# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the
# surface
VlnPlot(study2, "adt_CD19") #Cluster 5 is showing up! We could see this visually from previous plots

# we can also identify alternative protein and RNA markers for this cluster through differential
# expression
adt_markers <- FindMarkers(study2, ident.1 = 5, assay = "ADT")
rna_markers <- FindMarkers(study2, ident.1 = 5, assay = "RNA")

head(adt_markers)
head(rna_markers)




#------------------------------------------------------------------------#
#----------------STEP 6: Additional Visualizations-----------------------#
#------------------------------------------------------------------------#

# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if
# desired by using HoverLocator and FeatureLocator
FeatureScatter(study2, feature1 = "adt_CD19", feature2 = "adt_CD3")

# view relationship between protein and RNA
FeatureScatter(study2, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(study2, feature1 = "adt_CD4", feature2 = "adt_CD8")




#------------------------------------------------------------------------#
#----------------STEP 7: Save the Seurat Object--------------------------#
#------------------------------------------------------------------------#

saveRDS(study2, paste0(outputLoc, "SeuratObjects/seuratObj_ERP124005_subject1.rda"))


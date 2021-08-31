# MOSCATO - Leukemia Example

These files reproduce a real data example using MOSCATO (supporting code shared at https://github.com/lorinmil/MOSCATO). This real data example combines multi-omic single-cell data generated by CITE-seq technology for healthy patients and patients with leukemia. The outcome is leukemia versus healthy, and the two data types are RNA (i.e., X) and ADT (i.e., G). Before applying MOSCATO, the data must first be properly downloaded and pre-processed. The data must be integrated across subjects and cells clustered which is done using Seurat version 4.0.3 (proposed by Hao et al in 2021). There are a total of 21 subjects: 14 healthy and 7 with leukemia.

## Downloading and pre-processing the subject-level data

### Study 1

The first study used was ERP124005, and the data was downloaded at https://data.humancellatlas.org/explore/projects/efea6426-510a-4b60-9a19-277e52bfa815/project-matrices. Files Seurat_ERP124005_Sample1.R through Seurat_ERP124005_Sample10.R provide the R code to pre-process the subject-level data individually.

### Study 2

The second study used was GSE152469, and the data was downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4616298. Seurat_GSE152469.R provides the R code to pre-process the subject from this study.

### Study 3

The last study used was GSE139369, and the data was downloaded at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139369 where the samples used were GSM4138872, GSM4138874, GSM4138875, GSM4138876, GSM4138879, GSM4138880, GSM4138881, GSM4138883, GSM4138885, and GSM4138886. The files starting with Seurat_GSE139369_ provide the R code to pre-process the subject-level data individually.

## Integrate the data across subjects

Seurat_Combine_All.R provides the R code to integrate and cluster the cells across the subjects.


# TME analysis 
# this program is used to re-inte T cells and Myeloid cells
# and to define cell types 
# 2021-12-29
#===========================================================================================================================================
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")

GBM <- readRDS("")
#===========================================================================================================================================
# re-cluster for Myeloids

Myeloid <- subset(tmp.dat,cells=which(tmp.dat$celltype=="Myeloid"))
inte.list <- list()
samplelist <- unique(Myeloid$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(Myeloid,cells=which(Myeloid$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=50)
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte,resolution=1)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_Myeloid.RDS")































#========================================================================================================================================================
# for Lymphcyte 
# 2021-12-30
#========================================================================================================================================================
library(Seurat)















































































































































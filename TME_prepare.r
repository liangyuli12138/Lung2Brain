
# analysis TME and Myeloid cells
# 2021-5-20
###############################################################################################
# 0. some annoatation about inte TME # use merge is also OK,because we do not use this assay
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")
# GBM1 <- readRDS("/public/workspace/lily/PS/Final_716/lesion1.RDS")
# GBM2 <- readRDS("/public/workspace/lily/PS/Final_716/lesion2.RDS")
# # integration samples 
# integration.anchors <- FindIntegrationAnchors(object.list = c(dat,GBM1,GBM2))
# inte <- IntegrateData(anchorset = integration.anchors)
# #FindVariableFeatures
# inte <- FindVariableFeatures(inte)
# ##Scaling the integrateda
# all.genes <- rownames(inte)
# inte <- ScaleData(inte, features = all.genes)
# #PCA
# inte <- RunPCA(inte)
# #cluster
# inte <- FindNeighbors(inte)
# inte <- FindClusters(inte)
# #TSNE
# # if Umap can not use
# inte <- RunTSNE(inte)

# #=================================================================================================================================
# # add some information 
# #
# #=================================================================================================================================
# inte$type_group[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- "GBM" # sample info 
# # cell type info 
# inte$type[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- inte$celltype[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))]
# inte$type[which(inte$type%in%c("BMDM","MG"))] <- "Myeloid"
# inte$type[which(inte$type=="Tumor_cell")] <- "maliganant"

# saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")

# 1. subset non-tumor samples and re-inte 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
sub.dat <- subset(dat,cells=which(dat$malignant=="non-tumor"))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_ntumor.RDS")

# re - inte 
#===================================================================================================================================
# 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_ntumor.RDS")
inte.list <- list()
samplelist <- unique(dat$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
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
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_ntumor.RDS")

# subset T,B cell and Myeloid cells 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_ntumor.RDS")
sub.lym <- subset(dat,cells=which(dat$type%in%c("B_cell","T_cell")))
sub.mye <- subset(dat,cells=which(dat$type%in%c("Myeloid")))

saveRDS(sub.lym,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")

saveRDS(sub.mye,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_myeloid.RDS")




########################################################################################################################################
# 2021-6-2
# analysis lymphocyte 
#=======================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")
DefaultAssay(dat) <- "RNA"

inte.list <- list()
samplelist <- unique(dat$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=20)
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
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")



#=========================================================================================================================
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/T_cell_featureplot.pdf")
DefaultAssay(dat) <- "RNA"
DimPlot(dat,label=T)
DimPlot(dat,group.by="type_group")
FeaturePlot(dat,features=c("MS4A1","CD79A","CD79B"),label=T,order=F) # B cell 
FeaturePlot(dat,features=c("PTPRC","CD3D","CD3E"),label=T,order=F) # CD4 T
FeaturePlot(dat,features=c("CD8A","CD8B","CD3D"),label=T,order=F) # CD8 T 
FeaturePlot(dat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),label=T,order=F) # Treg 
FeaturePlot(dat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),label=T,order=F) # exhausted
FeaturePlot(dat,features=c("TCF7","SELL","LEF1","CCR7"),label=T,order=F) # naive T
FeaturePlot(dat,features=c("IRF4", "CREM", "NR4A2"),label=T,order=F) # Th17 
dev.off()


pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/T_cell_vlnplot.pdf",width=10)
DefaultAssay(dat) <- "RNA"
VlnPlot(dat,features=c("MS4A1","CD79A","CD79B"),pt.size=0) # B cell 
VlnPlot(dat,features=c("PTPRC","CD3D","CD3E"),pt.size=0) # CD4 T
VlnPlot(dat,features=c("GZMA","IL2","GZMK","GNLY","GZMB","IFNG"),pt.size=0) # cytotoxic T 
VlnPlot(dat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),pt.size=0) # Treg 
VlnPlot(dat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),pt.size=0) # exhausted
VlnPlot(dat,features=c("TCF7","SELL","LEF1","CCR7"),pt.size=0) # naive T
VlnPlot(dat,features=c("IRF4", "CREM", "NR4A2"),pt.size=0) # Th17 
VlnPlot(dat,features=c("MAF", "CXCR5", "CXCL13"),pt.size=0) # Th
dev.off()

# combine with singleR
library(Seurat)
library(SingleR)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")
ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")

# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_MonacoImmuneData.RDS")
# res.cluster<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


# 2021-6-7
#============================================================================================================
# define cell 
dat$hbmarker <- NULL
dat$llymarker <- NULL
dat$celltype <- NULL
dat$putativeTumor3 <- NULL
dat$res.2 <- NULL
dat$Phase <- NULL
dat$G2M.Score <- NULL
dat$S.Score <- NULL
dat$percent.mito <- NULL
dat$nUMI <- NULL
dat$nGene <- NULL
dat$RNA_snn_res.2 <- NULL

dat$celltype <- "Undefine"
dat$celltype[which(dat$seurat_clusters%in%c(12,14))] <- "Plasma" 
dat$celltype[which(dat$seurat_clusters%in%c(4))] <- "Bcell" #CD79A CD79B
dat$celltype[which(dat$seurat_clusters%in%c(11,16))] <- "NKcell" # KLRD 
dat$celltype[which(dat$seurat_clusters%in%c(6,7))] <- "Treg" # IL2RA FOXP3 
dat$celltype[which(dat$seurat_clusters%in%c(0,8,10))] <- "CD8 exhausted/cytotoxic" # PDCD1 LAG3 CTLA4
dat$celltype[which(dat$seurat_clusters%in%c(1,5))] <- "CD8 cytotoxic" # GZMA 
dat$celltype[which(dat$seurat_clusters%in%c(2,3))] <- "CD4 naive" # IL7R CCR7
dat$celltype[which(dat$seurat_clusters%in%c(13,15))] <- "Tfh" # CXCL13
dat$celltype[which(dat$seurat_clusters%in%c(9))] <- "Undefine"  # CD3D and CD3E no violin 

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_lympho.RDS")
# mv into Tcell file 







































# 2021-12-23
# prepare to RDS has OK in material/code_record.R
# this program is used to make RDS for each sample
# and then integration for each sample
# last for cell type annotation 
# 2022-2-14 remember : inteV6 is used Nature medicine data and Inte 16 is used NC data
#=============================================================================================================================================
# read all sample and normalized and integration

# library(Seurat)
# inte.list<- list()
# samplelist <- c("A20190305","A20190312","T-Bsc1_MT10","D0927","E0927","Pair_BM","Pair_LUNG","scrBT1431m","scrBT1432m","BT1296","BT1297","BT1291","BT1292")
# for(i in 1:length(samplelist)){
#     tmp_dat <- readRDS(paste0("/public/workspace/lily/Lung2Brain/RDS/",samplelist[i],".rds"))
#     tmp_dat = NormalizeData(object = tmp_dat)
# 	tmp_dat <- FindVariableFeatures(object = tmp_dat)
# 	# scaling
# 	all.genes <- rownames(x = tmp_dat)
# 	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
#     inte.list[[i]] <- tmp_dat
# }

# # now integration 
# integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
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
# inte <- RunUMAP(inte,dims=1:10)

# saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Data/inte_V6.RDS")










# this program is used to make RDS for each sample
# and then integration for each sample
# last for cell type annotation 
# 2021-12-30 
# test found SangSung data show better result 
#=============================================================================================================================================
# read all sample and normalized and integration
# scp from .3 to .5 then integration
# scp lily@202.195.187.3:/public/workspace/lily/Lung2Brain/HG38_Data/RDS/\{A20190305.RDS,A20190312.RDS,D0927.RDS,E0927.RDS,Pair_BM.RDS,Pair_LUNG.RDS,T-Bsc1_MT10.RDS\} ./

library(Seurat)
inte.list<- list()
samplelist <- grep("*.RDS$",dir("~/Lung2Brain/Version6/Prepare_Data"),value=T)
for(i in 1:length(samplelist)){
    tmp_dat <- readRDS(paste0("~/Lung2Brain/Version6/Prepare_Data/",samplelist[i]))
    # tmp_dat = NormalizeData(object = tmp_dat)
	# tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# # scaling
	# all.genes <- rownames(x = tmp_dat)
	# tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
    inte.list[[i]] <- tmp_dat
}

# now integration 
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
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")





















#============================================================================================================================================
# 2021-12-25 
# merry Christmas !
# get data from .3 and now define cell types
#============================================================================================================================================
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")

# DefaultAssay(dat) <- "RNA"
# pdf("/public/workspace/lily/Lung2Brain/Version6/Prepare/inteV6_marker_featureplot.pdf")
# DimPlot(dat,label=T,label.size=6)
# FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,label.size=3,order=T) # T cell
# FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),label=T,label.size=3,order=T) # Myeloid
# FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,label.size=3,order=T) # B cell 
# FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,label.size=3,order=T) # Oligodendrocyte
# FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,label.size=3,order=T) # Fibroblast/Vascular
# FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,label.size=3,order=T) # Endothelial
# FeaturePlot(dat,features=c("EGFR","EPCAM","KRAS"),label=T,label.size=3,order=T,cols=c("lightgrey", "red"))
# dev.off()

# pdf("/public/workspace/lily/Lung2Brain/Version6/Prepare/inteV6_marker_Vlnplot.pdf",width=20)
# DimPlot(dat,label=T,label.size=6)
# VlnPlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),pt.size=0) # T cell
# VlnPlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),pt.size=0) # Myeloid
# VlnPlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),pt.size=0) # B cell 
# VlnPlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),pt.size=0) # Oligodendrocyte
# VlnPlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),pt.size=0) # Fibroblast/Vascular
# VlnPlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),pt.size=0) # Endothelial
# VlnPlot(dat,features=c("EGFR","EPCAM","KRAS"),pt.size=0)
# dev.off()








#============================================================================================================================================
# 2021-12-30
# define cell type and prepare for infercnv
# 
#============================================================================================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

DefaultAssay(dat) <- "RNA"
pdf("/public/workspace/lily/Lung2Brain/Version6/Prepare/inteS16_marker_featureplot.pdf")
DimPlot(dat,label=T,label.size=6)
FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,label.size=3,order=T) # T cell
FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),label=T,label.size=3,order=T) # Myeloid
FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,label.size=3,order=T) # B cell 
FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,label.size=3,order=T) # Oligodendrocyte
FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,label.size=3,order=T) # Fibroblast/Vascular
FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,label.size=3,order=T) # Endothelial
FeaturePlot(dat,features=c("EGFR","EPCAM","KRAS"),label=T,label.size=3,order=T,cols=c("lightgrey", "red"))
dev.off()

pdf("/public/workspace/lily/Lung2Brain/Version6/Prepare/inteS16_marker_Vlnplot.pdf",width=20)
DimPlot(dat,label=T,label.size=6)
VlnPlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),pt.size=0) # T cell
VlnPlot(dat,features=c('CD19','CD68','FCGR3A','LYZ'),pt.size=0) # Myeloid
VlnPlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),pt.size=0) # B cell 
VlnPlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),pt.size=0) # Oligodendrocyte
VlnPlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),pt.size=0) # Fibroblast/Vascular
VlnPlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),pt.size=0) # Endothelial
VlnPlot(dat,features=c("EGFR","EPCAM","KRAS"),pt.size=0)
dev.off()








#============================================================================================================================================
# set celltype and sample group
# 2021-12-27
# C4 maybe normal cells and C9 maybe mixture cell clusters 
# C22 粒细胞 C24 Mast cells
#============================================================================================================================================
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")

# dat$celltype <- "unknow"
# dat$celltype[which(dat$seurat_clusters%in%c(0,1,21,23))] <- "Tcell"
# dat$celltype[which(dat$seurat_clusters%in%c(3,6,8,11,12,14,19,22,24))] <- "Myeloid"
# dat$celltype[which(dat$seurat_clusters%in%c(15,18))] <- "Bcell"
# dat$celltype[which(dat$seurat_clusters%in%c(7,10,16,9))] <- "Fibroblast"
# dat$celltype[which(dat$seurat_clusters%in%c(17))] <- "Oligo."
# dat$celltype[which(dat$seurat_clusters%in%c(20))] <- "Endothelial"
# dat$celltype[which(dat$seurat_clusters%in%c(2,5,13,4))] <- "Epithelial"


# # also add type_group information
# #============================================================================================================================================
# dat$type_group <- "unknow"
# dat$type_group[which(dat$orig.ident%in%c("A20190305","A20190312","D0927","T_Bsc1","E0927","Pair_BM"))] <- "LCBM"
# dat$type_group[which(dat$orig.ident%in%c("BT1291","BT1292","BT1296","BT1297","scrBT1431m","scrBT1432m"))] <- "nMLUAD"
# dat$type_group[which(dat$orig.ident%in%c("Pair_LUNG"))] <- "MLUAD"

# saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")












#============================================================================================================================================
# set celltype and sample group
# 2021-12-30
# C18 is Mast cells
#============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

dat$celltype <- "unknow"
dat$celltype[which(dat$seurat_clusters%in%c(3,4,8,15,24))] <- "Tcell"
dat$celltype[which(dat$seurat_clusters%in%c(0,5,7,9,11,13,19,21,22,25))] <- "Myeloid"
dat$celltype[which(dat$seurat_clusters%in%c(12,17))] <- "Bcell"
dat$celltype[which(dat$seurat_clusters%in%c(1,10,14,18,27))] <- "Fibroblast"
dat$celltype[which(dat$seurat_clusters%in%c(26))] <- "Oligo."
dat$celltype[which(dat$seurat_clusters%in%c(23))] <- "Endothelial"
dat$celltype[which(dat$seurat_clusters%in%c(2,6,16,20))] <- "Epithelial"


# also add type_group information
#============================================================================================================================================
dat$orig.ident[which(dat$orig.ident=="GSE131907")] <- dat$Sample[which(dat$orig.ident=="GSE131907")]
dat$type_group <- "unknow"
dat$type_group[which(dat$orig.ident%in%c("A20190305","A20190312","D0927","T_Bsc1","E0927","Pair_BM"))] <- "LCBM"
dat$type_group[which(dat$orig.ident%in%c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T25",'LUNG_T30','LUNG_T34'))] <- "nMLUAD"
dat$type_group[which(dat$orig.ident%in%c("Pair_Lung"))] <- "MLUAD"

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")





















































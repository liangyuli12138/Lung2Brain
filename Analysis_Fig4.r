
# 2022-2-18
# analysis about TME-Meylodid
# 1. add GBM data and just merge for all data
# 2. subset myeloid cells and T cells 
# 3. Identify Myeloid cell subsets 
#=================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
# mye <- subset(dat,cells=which(dat$celltype=="Myeloid"))
inte.list <- list() 
samplelist <- unique(dat$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
	inte.list[[i]] <- tmp
	names(inte.list)[i] <- samplelist[i]
}

# get GBM data 
gbm.sample <- grep("\\.RDS$",dir("~/Lung2Brain/Version6/GBM_Data"),ignore.case=T,value=T)
gbm.list <- list()
for(i in 1:length(gbm.sample)){
	tmp <- readRDS(paste0("~/Lung2Brain/Version6/GBM_Data/",gbm.sample[i]))
	DefaultAssay(tmp) <- "RNA"
	gbm.list[[i]] <- tmp
	names(gbm.list)[i] <- gbm.sample[i]
}

tmp.list <- c(inte.list,gbm.list)
tmp.dat <- tmp.list[[1]]
tmp.list[[1]] <- NULL

tmp.dat <- merge(tmp.dat,tmp.list)

# set cell type information

# 1. sample info 
tmp.dat$sample <- tmp.dat$orig.ident
# tmp.dat$sample[which(tmp.dat$sample=="P673_jc")] <- "P673"
# tmp.dat$sample[which(tmp.dat$sample=="P689_2")] <- "P689"
# tmp.dat$sample[which(tmp.dat$sample=="P912_left0514")] <- "P912_L"
# tmp.dat$sample[which(tmp.dat$sample=="P912_right0514")] <- "P912_R"
tmp.dat$sample[which(tmp.dat$sample=="Pair_BM")] <- "PLCBM1"
tmp.dat$sample[which(tmp.dat$sample=="Pair_Lung")] <- "PLUNG1"
tmp.dat$sample[which(tmp.dat$sample=="RD-20180817-001-SR18271")] <- "RD001"
tmp.dat$sample[which(tmp.dat$sample=="RD-20180817-002-SR18271")] <- "RD002"


# 2. type_group info 
tmp.dat$type_group[which(tmp.dat$sample%in%c("h3","h4","RD001","RD002"))] <- "GBM"


# 3. cell type info 
tmp.dat$celltype.refine[which(tmp.dat$sample=="h3")] <- tmp.dat$celltype4[which(tmp.dat$sample=="h3")]
tmp.dat$celltype.refine[which(tmp.dat$sample=="h4")] <- tmp.dat$celltype4[which(tmp.dat$sample=="h4")]
tmp.dat$celltype.refine[which(tmp.dat$sample=="RD001")] <- tmp.dat$llymarker[which(tmp.dat$sample=="RD001")]
tmp.dat$celltype.refine[which(tmp.dat$sample=="RD002")] <- tmp.dat$llymarker[which(tmp.dat$sample=="RD002")]
# tmp.dat$celltype.refine[which(tmp.dat$sample=="P673")] <- tmp.dat$CCC.type[which(tmp.dat$sample=="P673")]
# tmp.dat$celltype.refine[which(tmp.dat$sample=="P689")] <- tmp.dat$CCC.type[which(tmp.dat$sample=="P689")]
# tmp.dat$celltype.refine[which(tmp.dat$sample=="P912_L")] <- tmp.dat$CCC.type[which(tmp.dat$sample=="P912_L")]
# tmp.dat$celltype.refine[which(tmp.dat$sample=="P912_R")] <- tmp.dat$CCC.type[which(tmp.dat$sample=="P912_R")]


# do some change 
tmp.dat$celltype.refine[which(tmp.dat$celltype.refine%in%c("fibroblast_vascular","Fibroblast/Vascular"))] <- "Fibroblast"
tmp.dat$celltype.refine[which(tmp.dat$celltype.refine%in%c("myeloid","Macrophage"))] <- "Myeloid"
tmp.dat$celltype.refine[which(tmp.dat$celltype.refine%in%c("oligodendrocyte","Oligodendrocyte"))] <- "Oligo."
tmp.dat$celltype.refine[which(tmp.dat$celltype.refine%in%c("T_cell","Tcells"))] <- "Tcell"
tmp.dat$celltype.refine[which(tmp.dat$celltype.refine%in%c("Tumor cell","malignant"))] <- "Tumor"

#==================================================================================================================================
saveRDS(tmp.dat,file="")


# subset myeloid data and T cell data
# these codes are for analysis with myeloids
myeloid <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Myeloid"))
#===================================================================================================================================
#
# this code is used for integration data 
# inte.list <- list() 
# samplelist <- unique(myeloid$sample)
# for(i in 1:length(samplelist)){
# 	tmp <- subset(dat,cells=which(dat$sample==samplelist[i]))
#     DefaultAssay(tmp) <- "RNA"
# 	inte.list[[i]] <- tmp
# 	names(inte.list)[i] <- samplelist[i]
# }

library(harmony)
DefaultAssay(myeloid) <- "RNA"
myeloid <- ScaleData(object = myeloid)
myeloid <- FindVariableFeatures(object = myeloid)
myeloid <- RunPCA(object = myeloid, features = VariableFeatures(object = myeloid))
myeloid <- RunHarmony(myeloid,group.by.vars="sample")
myeloid <- RunUMAP(myeloid,reduction="harmony",dims=1:30)
myeloid <- FindNeighbors(myeloid,reduction="harmony",dims=1:30)
myeloid <- FindClusters(myeloid)


# define small celltype 
pdf("~/tmp.pdf")
FeaturePlot(myeloid,features=c("CTSS","FCN1","S100A8","S100A9","LYZ","VCAN"),label=T,label.size=3) # monocyte 
FeaturePlot(myeloid,features=c("LGMN","CTSB","CD14","FCGR3A"),label=T,label.size=3) # macrophage
FeaturePlot(myeloid,features=c("MARCO","FABP4","MCEMP1"),label=T,label.size=3) # alveolar macrophage 
FeaturePlot(myeloid,features=c("CLEC10A","CD1C","CLEC4C","PTCRA","CCR7","LAMP3"),label=T,label.size=3) # DC
dev.off()


pdf("~/tmp.pdf",width=15)
VlnPlot(myeloid,features=c("CTSS","FCN1","S100A8","S100A9","LYZ","VCAN"),pt.size=0) # monocyte 
VlnPlot(myeloid,features=c("LGMN","CTSB","CD14","FCGR3A"),pt.size=0) # macrophage
VlnPlot(myeloid,features=c("MARCO","FABP4","MCEMP1"),pt.size=0) # alveolar macrophage 
VlnPlot(myeloid,features=c("CLEC10A","CD1C","CLEC4C","PTCRA","CCR7","LAMP3"),pt.size=0) # DC
dev.off()












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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_Myeloid_16s.RDS")




























































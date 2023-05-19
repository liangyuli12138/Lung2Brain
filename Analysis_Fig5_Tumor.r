
# this program was designed to analysis tumor cells in LCBM
# 1. velocyto for all LCBM tumor cells


#====================================================================================================================================================
# 1. velocyto for all tumor cells 
# get tumor cell looms for each samples
bytlib load libraries/hdf5-1.8.13
bytlib load R-3.6.0

library(loomR)
library(pagoda2)
library(SCopeLoomR)
library(SeuratWrappers)
library(rlist)
library(Seurat)

tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
tmp.dat$cellname <- sapply(strsplit(colnames(tmp.dat),"_"),function(x){x[[1]]})

rm(list=ls()[-which(ls()=="tmp.dat")])
samplename = "T_Bsc1"
LCBM<- as.Seurat(ReadVelocity(paste0("~/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/",samplename,".loom")))
LCBM$cellname <- gsub("x$|^3cellranger:","",colnames(LCBM))
LCBM.obj <- subset(tmp.dat,cells=which(tmp.dat$orig.ident==samplename&tmp.dat$celltype.refine=="Tumor"))

# # pair LCBM samples
# samplename = "R21125541_hg38"
# LCBM<- as.Seurat(ReadVelocity(paste0("~/Lung2Brain/Version6/Data/Velocyto/",samplename,".loom")))
# LCBM$cellname <- gsub("x$|^R21125541_hg38:","",colnames(LCBM))
# LCBM.obj <- subset(tmp.dat,cells=which(tmp.dat$orig.ident=="Pair_BM"&tmp.dat$celltype.refine=="Tumor"))
# LCBM.subset <- subset(LCBM,cells=which(LCBM$cellname%in%LCBM.obj$cellname))
# ncol(LCBM.subset)==ncol(LCBM.obj)
# LCBM.subset$orig.ident <- "PairLCBM"
# LCBM.subset <- RenameCells(LCBM.subset,new.names=paste0(LCBM.subset$orig.ident,"_",LCBM.subset$cellname))
# saveRDS(LCBM.subset,file=paste0("~/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/PairLCBM_tumor.loom",".RDS"))


# get tumor cell loom file 
LCBM.subset <- subset(LCBM,cells=which(LCBM$cellname%in%LCBM.obj$cellname))
ncol(LCBM.subset)==ncol(LCBM.obj)
LCBM.subset$orig.ident <- samplename
LCBM.subset <- RenameCells(LCBM.subset,new.names=paste0(LCBM.subset$orig.ident,"_",LCBM.subset$cellname))
saveRDS(LCBM.subset,file=paste0("~/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/",samplename,"_tumor.loom",".RDS"))


# merge all tumor cells loom file 
pair_LCBM <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/PairLCBM_tumor.loom.RDS")
A05 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/A20190305_tumor.loom.RDS")
A12 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/A20190312_tumor.loom.RDS")
D0927 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/D0927_tumor.loom.RDS")
E0927 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/E0927_tumor.loom.RDS")
T_Bsc1 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/T_Bsc1_tumor.loom.RDS")


loomdat <- merge(x=pair_LCBM,y=c(A05,A12,D0927,E0927,T_Bsc1))
loomdat@active.assay <- "ambiguous"
saveRDS(loomdat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/Merge_tumor.loom.rds")


# Seurat V4 to trans
# r/4.1.2
library(Seurat)
library(SeuratDisk)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/Merge_tumor.loom.rds")
dat@active.assay <- "ambiguous"
SaveH5Seurat(dat, filename = "/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/Merge_tumor_loom.h5Seurat")
Convert("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/Merge_tumor_loom.h5Seurat", dest = "h5ad")




#======================================================================================================================================================
# get scvelo res and calculate
library(rhdf5)
h5f = H5Fopen("~/Lung2Brain/Version6/Data/Velocyto/Res/Merge_tumor_scvelo_res.h5ad")
tmp.dat <- data.frame(cellname=h5f$obs$`_index`,
  latent_time=h5f$obs$latent_time,pseudotime=h5f$obs$velocity_pseudotime,root_cell=h5f$obs$root_cells)

saveRDS(tmp.dat,file="~/Lung2Brain/Version6/Data/Velocyto/Res/Merge_tumor_scvelo_res_time.RDS")






# analysis BMS and scvelo latent time
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM"))

dat$cellname <- paste0(dat$orig.ident,"_",sapply(strsplit(colnames(dat),"_"),function(x){x[[1]]}))
dat$cellname <- gsub("^Pair_BM","PairLCBM",dat$cellname)
# length(which(dat$cellname%in%scvelo_res$cellname))
dat <- RenameCells(dat,new.names=dat$cellname)

scvelo_res <- readRDS("~/Lung2Brain/Version6/Data/Velocyto/Res/Merge_tumor_scvelo_res_time.RDS")
rownames(scvelo_res) <- scvelo_res$cellname
scvelo_res <- scvelo_res[colnames(dat),]
all(colnames(dat)==rownames(scvelo_res))



# 
gene <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Signature/BMS_V6_gene.RDS")
dat <- AddModuleScore(dat,features=list(gene),name="BMS")
cor.test(dat$BMS1,scvelo_res$latent_time)
cor.test(dat$BMS1,scvelo_res$latent_time,method="spearman")


# ssGSEA calculate
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_V6"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
cor.test(mod[,2],scvelo_res$pseudotime,method="spearman")
dat$BMS <- mod[,2]


# set BMS group 
dat$BMS_group <- "Medium"
dat$BMS_group[which(dat$BMS>quantile(dat$BMS,0.67))] <- "BMSH"
dat$BMS_group[which(dat$BMS<quantile(dat$BMS,0.33))] <- "BMSL"
library(future)
future::plan(multisession, workers=20)
markerl <- FindMarkers(dat,assay="RNA",group.by="BMS_group",ident.1="BMSL",ident.2="BMSH",logfc.threshold=0,min.pct=0.05)













#=================================================================================================================================================
# 2022-6-22
# re-integration and re-cluster all LCBM tumor cells
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM"))


options(future.globals.maxSize= 8912896000) # 8500MB change size
inte.list <- list()
samplelist <- unique(dat$orig.ident)
for(i in 1:length(samplelist)){
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
inte <- FindClusters(inte,resolution=1)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte.RDS")














#===========================================================================================================================
# 2022-6-29
# analysis CNMF

useage <- read.csv("/public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_cNMF/LCBM_cNMF.usages.k_4.dt_0_2.consensus.txt",sep="\t",stringsAsFactors=F)
rownames(useage) <- useage$X
useage$X <- NULL

genescore <- read.csv("/public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_cNMF/LCBM_cNMF.gene_spectra_score.k_4.dt_0_2.txt",sep="\t",stringsAsFactors=F)
genescore$X <- NULL
genescore <- t(genescore)
colnames(genescore) <- c("P1","P2","P3","P4")




library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte.RDS")

# program 
useage <- read.csv("/public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_cNMF/LCBM_cNMF.usages.k_4.dt_0_2.consensus.txt",sep="\t",stringsAsFactors=F)
rownames(useage) <- useage$X
useage$X <- NULL
pdat <- apply(useage,2,function(x){(x-min(x))/(max(x)-min(x))})

# change names
dat$cells <- gsub("_BM","LCBM",paste0(dat$orig.ident,"_",sapply(strsplit(colnames(dat),"_"),function(x){x[[1]]})))
length(which(dat$cells%in%rownames(pdat)))
pdat <- pdat[dat$cells,]
dat$P1 <- pdat[,1]
dat$P2 <- pdat[,2]
dat$P3 <- pdat[,3]
dat$P4 <- pdat[,4]
FeaturePlot(dat,features=c("P1","P2","P3","P4"))
FeaturePlot(dat,features=c("P1","P2","P3","P4"),order=T)

# calculate BMS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_V6","HPSC_C5"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte_BMS_mod.RDS")









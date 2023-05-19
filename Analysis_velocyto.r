# this program is used to analysis velocyto in R 


# use velocyte to get loom files
# run in 202.195.187.3
# Lung sample 
conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/GRCh38_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/lily/LCBM_pair/R21136163_hg38/ \
    /public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf


# velocyto run10x -m /public/workspace/lily/REF/GRCh38_rmsk.gtf \
#     -@ 10 --samtools-memory 2000 \
#     /public/workspace/zhumy/test/panyuqin_fibroblast/result/cellranger_res1/10XT2H2/ \
#     /public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf

# LCBM sample 
# need to cp  to /public/workspace/lily/LCBM_pair/R21125541_hg38
cp -r /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/WSY/3cellranger/* /public/workspace/lily/LCBM_pair/R21125541_hg38/
# 
conda activate velocity
module load  samtools-1.9

velocyto run10x -m /public/workspace/lily/REF/GRCh38_rmsk.gtf \
    -@ 10 --samtools-memory 2000 \
    /public/workspace/lily/LCBM_pair/R21125541_hg38/ \
    /public/workspace/xiapeng/CellRanger/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf



##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
# now run in R 



bytlib load libraries/hdf5-1.8.13
bytlib load R-3.6.0
R
#source('/public/workspace/caiyun/R/functions.R')
#source('/public/workspace/caiyun/R/requiredPcg.R')
library(loomR)
library(velocyto.R)
library(pagoda2)
library(SCopeLoomR)
library(SeuratWrappers)
library(rlist)
library(Seurat)

tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
LCBM<- as.Seurat(ReadVelocity("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/R21125541_hg38.loom"))
LCBM$cellname <- gsub("x$|^R21125541_hg38:","",colnames(LCBM))
LCBM.obj <- subset(tmp.dat,cells=which(tmp.dat$orig.ident=="Pair_BM"))
LCBM.obj$cellname <- gsub("\\_14$","",colnames(LCBM.obj))
tumor.LCBM <- subset(LCBM.obj,cells=which(LCBM.obj$celltype.refine=="Tumor"))

# subset cell 
LCBM.subset <- subset(LCBM,cells=which(LCBM$cellname%in%tumor.LCBM$cellname))
LCBM.subset = NormalizeData(LCBM.subset, verbose = FALSE)
LCBM.subset = FindVariableFeatures(LCBM.subset, selection.method = "vst", verbose = FALSE)
LCBM.subset <- RenameCells(LCBM.subset,new.names= paste0("LCBM_",LCBM.subset$cellname))

# some tumor cell do not found ,try to subset seurat obj 
tumor.LCBM.subset <- subset(tumor.LCBM,cells=which(tumor.LCBM$cellname%in%LCBM.subset$cellname))
tumor.LCBM.subset <- RenameCells(tumor.LCBM.subset,new.names=paste0("LCBM_",tumor.LCBM.subset$cellname))
saveRDS(tumor.LCBM.subset,file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_LCBM.tumor.RDS")



#==========================================================================================================================================================
# do the same prepare for paired Lung sample
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
Lung<- as.Seurat(ReadVelocity("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/R21136163_hg38.loom"))
Lung$cellname <- gsub("x$|^R21136163_hg38:","",colnames(Lung))
Lung.obj <- subset(tmp.dat,cells=which(tmp.dat$orig.ident=="Pair_Lung"))
Lung.obj$cellname <- gsub("\\_15$","",colnames(Lung.obj))
tumor.Lung <- subset(Lung.obj,cells=which(Lung.obj$celltype.refine=="Tumor"))

# subset cell 
Lung.subset <- subset(Lung,cells=which(Lung$cellname%in%tumor.Lung$cellname))
Lung.subset = NormalizeData(Lung.subset, verbose = FALSE)
Lung.subset = FindVariableFeatures(Lung.subset, selection.method = "vst", verbose = FALSE)
Lung.subset <- RenameCells(Lung.subset,new.names= paste0("Lung_",Lung.subset$cellname))

# some tumor cell do not found ,try to subset seurat obj 
tumor.Lung.subset <- subset(tumor.Lung,cells=which(tumor.Lung$cellname%in%Lung.subset$cellname))
tumor.Lung.subset <- RenameCells(tumor.Lung.subset,new.names=paste0("Lung_",tumor.Lung.subset$cellname))
saveRDS(tumor.Lung.subset,file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_Lung.tumor.RDS")


#==========================================================================================================================================================

loomDat <- merge(LCBM.subset,Lung.subset)
saveRDS(loomDat,file='/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.rds')

#==========================================================================================================================================================
# integration 
integration.anchors <- FindIntegrationAnchors(object.list = c(tumor.LCBM.subset,tumor.Lung.subset))
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
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_inte.tumor.RDS")


# merge result is not OK
#==========================================================================================================================================================
#==========================================================================================================================================================
#==========================================================================================================================================================
#==========================================================================================================================================================
#==========================================================================================================================================================
#==========================================================================================================================================================
#==========================================================================================================================================================
library(loomR)
library(velocyto.R)
library(pagoda2)
library(SCopeLoomR)
library(SeuratWrappers)
library(rlist)
library(Seurat)

savepath="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_res/"
dat_rdsPath ="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_inte.tumor.RDS"
loomDat_rdsPath ="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.rds"
spliced_minAvg=0.5
unspliced_minAvg=0.1
kCells=500
n=5000
fit_quantile=0.05
param="orig.ident"

get_color_scheme = function(type = "clusters") {
  library(ggsci)
  if (type == "samples") {
    color_scheme = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
  }
  if (type == "clusters") {
    color_scheme = c( pal_d3("category20")(20), pal_d3("category20b")(20), pal_d3("category20c")(20),pal_igv("default")(51))
  }
  return(color_scheme)
}

### mainText
# savepath = paste0(gsub('/$','',savepath),'/')
# if(!file.exists(savepath)) dir.create(savepath)

DAT = readRDS(dat_rdsPath)
loomDat = readRDS(loomDat_rdsPath)

DAT <- subset(DAT,cells=which(DAT$seurat_clusters%in%c(8,9,2,6,7)))
loomDat <- subset(loomDat,cells=colnames(DAT))

loomDat@assays$integrated = DAT@assays$integrated
loomDat@reductions = DAT@reductions
loomDat = AddMetaData(loomDat,DAT@meta.data)
loomDat$param = as.character(loomDat[[param]][,1])

cluster.colors = get_color_scheme('clusters')[1:length(unique(loomDat$param))]
names(cluster.colors) = sort(as.character(unique(loomDat$param)))
cell.colors = sapply(loomDat$param,function(x) {cluster.colors[as.character(x)]})
names(cell.colors) = colnames(loomDat)

loomDat = RunVelocity(object = loomDat, deltaT = 1, kCells = as.numeric(kCells), fit.quantile = as.numeric(fit_quantile),
                              spliced.average = as.numeric(spliced_minAvg),unspliced.average=as.numeric(unspliced_minAvg),
                              reduction = 'umap',group.by='param',ncores=8)
saveRDS(loomDat,file=paste0(savepath,param,"_",'loomDat.rds'))

library(ggplot2)
pdf(paste0(savepath,"cell_velocity_umap_",kCells,'kCells_',n,"_",param,".pdf"),useDingbats=F,width=10,height=10)
DimPlot(loomDat,reduction='umap',group.by='param',cols=cluster.colors) + theme(aspect.ratio=1)
dev.off()

pdf(paste0(savepath,"cell_velocity_",kCells,'kCells_',n,"_",param,".pdf"),useDingbats=F,width=10,height=10)
rs = show.velocity.on.embedding.cor(emb = Embeddings(object = loomDat, reduction = "umap"), 
     vel = Tool(object = loomDat, slot = "RunVelocity"), n = as.numeric(n), scale = "sqrt", 
     cell.colors = ac(x = cell.colors, alpha = 0.5), cex = 0.8, arrow.scale = 3, 
     show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
    do.par = FALSE, cell.border.alpha = 0.1,return.details=TRUE)
dev.off()

save = list(transitionProbability = rs$tp,arrowEstimatesPos = rs$arrows, scale = rs$scale)
saveRDS(save,file=paste0(savepath,param,"_",'velocyto.rds'))










# 2022-2-11
#==================================================================================================================================
# calculate RNA velocity 
savepath="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_res/"
dat_rdsPath ="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_inte.tumor.RDS"
loomDat_rdsPath ="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.rds"
spliced_minAvg=0.5
unspliced_minAvg=0.1
kCells=500
n=5000
fit_quantile=0.05
param="orig.ident"


# need to add loomDat
loomDat <- readRDS(paste0(savepath,param,"_",'loomDat.rds'))
DimPlot(loomDat,reduction='umap',group.by='seurat_clusters',label=T,label.size=10) + theme(aspect.ratio=1)

save <- readRDS(paste0(savepath,param,"_",'velocyto.rds'))
# dat <- save$arrowEstimatesPos[which(rownames(save$arrowEstimatesPos)%in%colnames(loomDat)[which(loomDat$seurat_clusters%in%c(0:9))]),]
dat <- save$arrowEstimatesPos

dat1 <- dat[grep("^LCBM",rownames(dat)),]
dat2 <- dat[grep("^Lung",rownames(dat)),]

Velocyto_cal <- function(dat1,dat2){
  # just calculate ponits in dat1 ,if want to calcualte dat2 ,change!
  # dat2 is just as ref data
  mean.point1 <- c(mean(dat1[,1]),mean(dat1[,2])) # use all start point 
  mean.point2 <- c(mean(dat2[,1]),mean(dat2[,2]))
  tmp.res <- data.frame(t(apply(dat1,1,function(x){
    Bdist1 <- sqrt((x[1]-mean.point1[1])^2 + (x[2]-mean.point1[2])^2) # before is use start - mean point (lesion1)
    Adist2 <- sqrt((x[3]-mean.point1[1])^2 + (x[4]-mean.point1[2])^2) # after is use end - mean point (lesion1)
    Bdist3 <- sqrt((x[1]-mean.point2[1])^2 + (x[2]-mean.point2[2])^2) # before is use start - mean point (lesion2/refernce data)
    Adist4 <- sqrt((x[3]-mean.point2[1])^2 + (x[4]-mean.point2[2])^2) # after is use end- mean point (lesion2/refernce data)
    dist.diff1 <- Bdist1 - Adist2
    dist.diff2 <- Bdist3 - Adist4

    c(Bdist1,Adist2,Bdist3,Adist4,dist.diff1,dist.diff2)
    })))
  colnames(tmp.res) <- c("Bdist1","Adist2","Bdist3","Adist4","dist.diff1","dist.diff2")

  return(tmp.res)
}

LCBM.d <- Velocyto_cal(dat1,dat2) # when calculate P673 ,should use P689 as ref
Lung.d <- Velocyto_cal(dat2,dat1)


LCBM.d.f <- LCBM.d[which((LCBM.d$dist.diff1 * LCBM.d$dist.diff2)<0),]
Lung.d.f <- Lung.d[which((Lung.d$dist.diff1 * Lung.d$dist.diff2)<0),]

# summary 
length(which(LCBM.d.f$dist.diff1<0&LCBM.d.f$dist.diff2>0)) # far LCBM , near Lung
length(which(LCBM.d.f$dist.diff1>0&LCBM.d.f$dist.diff2<0))
length(which(Lung.d.f$dist.diff1<0&Lung.d.f$dist.diff2>0))
length(which(Lung.d.f$dist.diff1>0&Lung.d.f$dist.diff2<0))

# 若方向是LCBM to Lung , 那么应该LCBM 是diff1<0,diff2>0 多； 而Lung应该是diff1>0,diff2<0 点更多
# OR = 11.6367 ; p-value=0.057







# 2022-6-8
# calculate velocity score 
# use this trans to hd5
# r/4.1.2
#====================================================================================================================
library(Seurat)
library(SeuratDisk)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.rds")
dat@active.assay <- "ambiguous"
SaveH5Seurat(dat, filename = "/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.h5Seurat")
Convert("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.h5Seurat", dest = "h5ad")



# Convert("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res.h5ad", dest = "h5seurat", overwrite = TRUE)
# tmp.dat <- LoadH5Seurat("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res.h5seurat")

# use hdf5 file to get info
# run in R-3.6.0
# https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
library(rhdf5)
h5f = H5Fopen("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res.h5ad")
tmp.dat <- data.frame(cellname=h5f$obs$cellname,
  latent_time=h5f$obs$latent_time,pseudotime=h5f$obs$velocity_pseudotime,root_cell=h5f$obs$root_cells)
saveRDS(tmp.dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res_time.RDS")












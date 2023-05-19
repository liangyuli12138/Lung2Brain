
#==============================lly new find markers

# 2020-6-5 lly used raw lung data to anlysis 
#======================================
#  functions


#=================================================================================================================
# BT1292 cellranger result is not good ,so lly have to use cellranger to rr-run 

#/public/workspace/zhumy/lung_brain/rawdata/E-MTAB-6149/BT1292


module load cellranger-3.0.2

cd /public/workspace/lily/wulx

cellranger count --id=T-Bsc1 \
   --fastqs=/public/workspace/wulx/ANNUO/data/ANCJS190333_PM-JS190333-01_AHYJWHCCXY_2019-02-26/Cleandata/T-Bsc1/ \
   --transcriptome=/public/workspace/zhumy/ref/refdata-cellranger-hg19-1.2.0 \
   --chemistry=SC3Pv3 \
   --localcores=20 --localmem=60 --nosecondary
   --expect-cells=10000


#==========Lung data
#======================================================================================= lung data 
# BT1292 BT1293 BT1297 BT1294 BT1432m BT1429m


name="scrBT1430m"
# scrBT1429m
# scrBT1432m
# scrBT1430m
datafile = paste0("/public/workspace/lily/Lung2Brain/Lung_Raw/E-MTAB-6653/",name,"/outs/filtered_feature_bc_matrix")
outpath = "/public/workspace/lily/Lung2Brain/Prepare/Lung/"
rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

tmp <- Read10X(data.dir = datafile)
dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#=======================filter cells 
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
saveRDS(dat,file=paste0(rdspath,name,'.rds'))


#============================== cell ranger result in zhumengy's file
# /public/workspace/zhumy/lung_brain/result/cellranger_result/E-MTAB-6149

sample <- c("BT1296","BT1297","BT1294","BT1295")
sample <- c("BT1299","BT1300")
for(i in 1:length(sample)){

	name <- "BT1299"
	datafile = paste0("/public/workspace/zhumy/lung_brain/result/cellranger_result/E-MTAB-6149/",name,"/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/Prepare/Lung/"
	rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

	tmp <- Read10X(data.dir = datafile)
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()



	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
	saveRDS(dat,file=paste0(rdspath,name,'.rds'))


}




#==========Brain Metastasis data
# Brain metastasis all use same subset (200,7500) and percent.mt <20
# outpath='/public/workspace/lily/Lung2Brain/Prepare/BM/'
name = 'T-Bsc1'
# A20190312
# A20190305

datafile = paste0("/public/workspace/zhumy/lung_brain/result/cellranger_result/",name,"/outs/filtered_feature_bc_matrix")
outpath = "/public/workspace/lily/Lung2Brain/Prepare/BM/"
rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

tmp <- Read10X(data.dir = datafile)
dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#==================need change in different samples 
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
saveRDS(dat,file=paste0(rdspath,name,"_MT10",'.rds'))



#==========GBM data
name = 'RD-20180817-002-SR18271'
# RD-20180817-001-SR18271
# 
# RD-20180820-001-SR18271


datafile = paste0("/public/workspace/wulx/SHIHE/result/",name,"/outs/filtered_gene_bc_matrices/hg19/")
outpath = "/public/workspace/lily/Lung2Brain/Prepare/GBM/"
rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

tmp <- Read10X(data.dir = datafile)
dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#==================need change in different samples 
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(dat,file=paste0(rdspath,name,'.rds'))






#==========================================================================================================================================================================================
#=====================================
# when all samples have filter cells ,you have to run PCA and find clusters

filepath <- "/public/workspace/lily/Lung2Brain/RDS/"
sample <- dir(filepath) #14

for(i in 1:length(sample)){
	dat <- readRDS(paste0(filepath,sample[i]))

# up step not normalize and find cluster
	dat = NormalizeData(object = dat)
	dat <- FindVariableFeatures(object = dat)
	all.genes <- rownames(x = dat)
	dat <- ScaleData(object = dat, features = all.genes)

	# Run PCA
	dat <- RunPCA(object = dat, features = VariableFeatures(object = dat),npcs = 50)
	pdf(paste0("/public/workspace/lily/Lung2Brain/Prepare/ElbowPlot_",gsub("\\.rds$","",sample[i]),".pdf"))
	p<- ElbowPlot(object = dat,ndims=50)
	print(p)
	dev.off()

	pdf(paste0("/public/workspace/lily/Lung2Brain/Prepare/DoHeatmap/",gsub("\\.rds$","",sample[i]),".pdf"))
	DimHeatmap(dat, dims = 1:10, cells = 500, balanced = TRUE)
	DimHeatmap(dat, dims = 11:20, cells = 500, balanced = TRUE)
	DimHeatmap(dat, dims = 21:30, cells = 500, balanced = TRUE)
	dev.off()


}

# by up steps ,lly think maybe PC=20 is the best choices 
# by use following code ,lly think maybe PC=10 is ok 



#=============================================================================
# findCluster 
# samples which cells numbers is smaller than 1W use resolution=1;which cells number is bigger than 1W use resolution=2

library(Seurat)
res1 = 1
res2 = 2

filepath <- "/public/workspace/lily/Lung2Brain/RDS/"
#sample <- dir(filepath) #14
sample <- c("BT1299.rds","BT1300.rds")
for(i in 1:length(sample)){
	# read data
	dat <- readRDS(paste0(filepath,sample[i]))

	# chose different resolution
	if(ncol(dat)>10000){ 
		tmp_res <- res2 
		}else{
			tmp_res <- res1
		}

# up step not normalize and find cluster
	dat = NormalizeData(object = dat)
	dat <- FindVariableFeatures(object = dat)
	all.genes <- rownames(x = dat)
	dat <- ScaleData(object = dat, features = all.genes)

	# Run PCA
	dat <- RunPCA(object = dat, features = VariableFeatures(object = dat),npcs = 50)

	# just caculate 50 PCs,but not all use 
	# lly thinks use 10 PCs is ok ,because as the Doheatmap function plot ,some PCs is not very ok ,and the tutor is use 10 .
	dat <- FindNeighbors(dat, dims = 1:20)
	dat <- FindClusters(dat, resolution = tmp_res)
	
	# the tutorial says should use the same dims in Findcluster function
	dat <- RunTSNE(object = dat,dims = 1:20)

	saveRDS(dat,file=paste0("/public/workspace/lily/Lung2Brain/Data/cluster_",sample[i]))


}



#===========================================================================================================================================================================================
# tumor identify 
library(Seurat)
library(pheatmap)
library(scCNA)
data(CpRMAP)

filepath <- "/public/workspace/lily/Lung2Brain/Data/"
sample <- dir(filepath)[grep("\\.rds$",dir(filepath))] #14

filepath <- "//public/workspace/lily/Lung2Brain/Data/cluster_data/"
sample <- c("cluster_BT1299.rds","cluster_BT1300.rds")

for(i in 1:length(sample)){
	# read data
	dat <- readRDS(paste0(filepath,sample[i]))

	dat@active.ident -> tag
	#names(tag) = rownames(info.f)
	dat[['RNA']]@data -> dat.f

	rs = identificationTumorCells(mat.f=log2(dat.f+1), cluster.tag=tag, chrpq=GRCh37.chrpq)

	#get the result information 
	save(rs,file=paste0('/public/workspace/lily/Lung2Brain/Data/rs/',gsub("\\.rds$","",sample[i]),"_rs.RData"))
	#rs$information$putativeTumor2
	dat <- AddMetaData(object = dat, metadata = rs$information$putativeTumor2, col.name = "putativeTumor2")

	dat <- AddMetaData(object = dat, metadata = rs$information$putativeTumor, col.name = "putativeTumor")

	saveRDS(dat, file = paste0('/public/workspace/lily/Lung2Brain/Data/tumoridentify/',gsub("\\.rds$","",sample[i]),"_identify.rds"))

}





#========================================================================================================================================================================
# inte all samples

library(Seurat)

inte.list <- list()
names <- c()
sample <- dir("/public/workspace/lily/Lung2Brain/Data/tumoridentify/")
for(i in 1:length(sample)){
	name <- gsub("^cluster_|.rds$","",sample[i])
	inte.list[[i]] <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Data/tumoridentify/",sample[i]))
	names <- c(names,name)
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list, dims = 1:30)
inte <- IntegrateData(anchorset = integration.anchors, dims = 1:30)

DefaultAssay(inte) <- "integrated"
#FindVariableFeatures
inte <- FindVariableFeatures(inte, selection.method = "vst")

##Scaling the integratedata
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)

saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte_rough.rds")
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte, dims = 1:20)
inte <- FindClusters(inte, resolution = 2)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte,dims = 1:20)
#save data
saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte.rds")





#================add sample info 
a = inte$orig.ident
#get 3 sample
a[which(a[]=="RD-20180817-001-SR18271")] <- 'S1701'
a[which(a[]=="RD-20180817-002-SR18271")] <- 'S1702'
a[which(a[]=="RD-20180820-001-SR18271")] <- 'S2001'

inte <- AddMetaData(object = inte, metadata = a, col.name = "sample_group")

#================correct tumor-nontumor cell 

inte$malignant <- inte$putativeTumor2
inte$malignant[which(inte$orig.ident=="BT1294")] <- "non-tumor"
inte$malignant[which(inte$orig.ident=="scrBT1429m")] <- "non-tumor"


# add group info 
dat$type_group <- "LUAD"
dat$type_group[which(dat$sample_group%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
dat$type_group[which(dat$sample_group%in%c("S1701","S1702","S2001"))]<- "GBM"


saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte.rds")



#=====================================================================
# get tumor cell and non-tumor cell
tumor <- subset(dat,cells=which(dat$malignant=="tumor"))
saveRDS(tumor,file="/public/workspace/lily/Lung2Brain/Data/run_data/tumor.rds")
ntumor <- subset(dat,cells=which(dat$malignant=="non-tumor"))
saveRDS(tumor,file="/public/workspace/lily/Lung2Brain/Data/run_data/non-tumor.rds")



#====================================================================
# tumor need to re-cluster
 tumor <- readRDS("/public/workspace/lily/Lung2Brain/Data/run_data/tumor.rds")
DefaultAssay(tumor) <- "integrated"

tumor <- FindVariableFeatures(object = tumor)
# scaling
all.genes <- rownames(x = tumor)
tumor <- ScaleData(object = tumor, features = all.genes)
# PCA
tumor <- RunPCA(object = tumor, features = VariableFeatures(object = tumor))
# clustering
tumor <- FindNeighbors(object = tumor,dims=1:20)
# select proper resolution
tumor <- FindClusters(object = tumor,resolution=3)
# T-SNE
tumor <- RunTSNE(object = tumor,dims=1:20)


#============================
# add group info 
tumor$type_group <- "LUAD"
tumor$type_group[which(tumor$sample_group%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
tumor$type_group[which(tumor$sample_group%in%c("S1701","S1702","S2001"))]<- "GBM"




pdf("./tmp.pdf")
DimPlot(tumor)
DimPlot(tumor,group.by="type_group",cols=c("#F19143","#5386E4","#4D8B31"))
DimPlot(tumor,label=T)
dev.off()






#=============================================================
# do not use edge tumor sample 
library(Seurat)

inte.list <- list()
names <- c()
sample <- dir("/public/workspace/lily/Lung2Brain/Data/tumoridentify/")
sample[which(sample[]=="cluster_BT1295_identify.rds")] <- NA
sample[which(sample[]=="cluster_scrBT1430m_identify.rds")] <- NA
sample <- na.omit(sample)
for(i in 1:length(sample)){
	name <- gsub("^cluster_|_identify\\.rds$","",sample[i])
	if(name!="BT1295"&name!="scrBT1430m"){
		inte.list[[i]] <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Data/tumoridentify/",sample[i]))
		names <- c(names,name)
	}
	
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list, dims = 1:30)
inte <- IntegrateData(anchorset = integration.anchors, dims = 1:30)

DefaultAssay(inte) <- "integrated"
#FindVariableFeatures
inte <- FindVariableFeatures(inte, selection.method = "vst")

##Scaling the integratedata
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)

saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte12_rough.rds")
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte, dims = 1:20)
inte <- FindClusters(inte, resolution = 2)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte,dims = 1:20)
#save data
saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte12.rds")

































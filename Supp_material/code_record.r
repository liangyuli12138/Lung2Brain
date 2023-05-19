

# code record 
# this code was renew in 2021-12-23
#==========================================================================================================================================
# 1. run cellranger 
# cell ranger was run by hg19,use cellranger v3.0.2
# cell ranger was run by zhumy and lly 
# lung cancer brain metastasis was run by zhumy
# sample TB sample A05 and sample A12
# new LCBM samples was run by lily 
# sample D0927 E0927 and paired LCBM 
# some lung cancer sample was run by lly
# sample scr1432m and sample scr1431m 
# and other lung cancer sample was run by zhumy 
# sample BT1291,BT1292,BT1296,BT1297

# 2. sample filter 
# Lung cancer data were filter by MT<20%, min.cells=3, min.features=200, nfeature_RNA>200 and nfeature_RNA<5000
############################################################################################################################################
############################################################################################################################################
#===========================================================================================================================================
# @ this sample do not need to normalized 
# cell ranger result in lily's file
# /public/workspace/lily/Lung2Brain/Lung_Raw/E-MTAB-6653/

sample <- c("scrBT1431m","scrBT1432m")
for(i in 1:length(sample)){

	name <- sample[i]
	datafile = paste0("/public/workspace/lily/Lung2Brain/Lung_Raw/E-MTAB-6653/",name,"/outs/filtered_feature_bc_matrix")
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

# cell ranger result in zhumengy's file
# /public/workspace/zhumy/lung_brain/result/cellranger_result/E-MTAB-6149

sample <- c("BT1296","BT1297","BT1291","BT1292")
for(i in 1:length(sample)){

	name <- sample[i]
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

#===========================================================================================================================================
# lung cancer brain metastasis data  # also run in 202.195.187.3
# Brain metastasis all use same subset (200,7500) and percent.mt <20
# ############################however TBSC use MT<10 because it has too many cells
# outpath='/public/workspace/lily/Lung2Brain/Prepare/BM/'
name = 'A20190305'
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

# maybe need to change in different samples #this is for TBSC
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
saveRDS(dat,file=paste0(rdspath,name,'.rds'))



#====================================================================================================================================
# new Brain Metastasis samples and Paired LCBM samples

sample <- c("D0927","E0927")
for(i in 1:length(sample)){

	name <- sample[i]
	datafile = paste0("/public/workspace/lily/Mutiple_LB/",name,"/",name,"/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/Prepare/BM/"
	rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

	tmp <- Read10X(data.dir = datafile)
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()

	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
	saveRDS(dat,file=paste0(rdspath,name,'.rds'))
}


# for LCBM pair 
sample <- c("R21125541","R21136163")
# R21136163 is LUNG 
# R21125541 is BM
for(i in 1:length(sample)){

	name <- sample[i]
	datafile = paste0("/public/workspace/lily/LCBM_pair/",name,"/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/Prepare/BM/"
	rdspath = '/public/workspace/lily/Lung2Brain/RDS/'

	if(sample[i]=="R21125541"){
		name <- "Pair_BM"
	}
	if(sample[i]=="R21136163"){
		name <- "Pair_LUNG"
	}
	tmp <- Read10X(data.dir = datafile)
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()

	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
	saveRDS(dat,file=paste0(rdspath,name,'.rds'))
}


















#=============================================================================================================================================
# GBM data was use PS data ,which would be used in TME analysis  # also run in 
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

# maybe need to change in different samples 
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
saveRDS(dat,file=paste0(rdspath,name,'.rds'))




#############################################################################################################################################
# Find cluster  and use scan to identify tumor cell 
#############################################################################################################################################
# findCluster 
# samples which cells numbers is smaller than 1W use resolution=1;which cells number is bigger than 1W use resolution=2

# library(Seurat)
# res1 = 1
# res2 = 2

# filepath <- "/public/workspace/lily/Lung2Brain/RDS/"
# #sample <- dir(filepath) #14
# sample <- c("BT1299.rds","BT1300.rds")
# for(i in 1:length(sample)){
# 	# read data
# 	dat <- readRDS(paste0(filepath,sample[i]))

# 	# chose different resolution
# 	if(ncol(dat)>10000){ 
# 		tmp_res <- res2 
# 		}else{
# 			tmp_res <- res1
# 		}

# # up step not normalize and find cluster
# 	dat = NormalizeData(object = dat)
# 	dat <- FindVariableFeatures(object = dat)
# 	all.genes <- rownames(x = dat)
# 	dat <- ScaleData(object = dat, features = all.genes)

# 	# Run PCA
# 	dat <- RunPCA(object = dat, features = VariableFeatures(object = dat),npcs = 50)

# 	# just caculate 50 PCs,but not all use 
#     # call.string: chr "RunTSNE(object = dat, dims = 1:20)
# 	dat <- FindNeighbors(dat, dims = 1:20)
# 	dat <- FindClusters(dat, resolution = tmp_res)
	
# 	# the tutorial says should use the same dims in Findcluster function
# 	dat <- RunTSNE(object = dat,dims = 1:20)

# 	saveRDS(dat,file=paste0("/public/workspace/lily/Lung2Brain/Data/cluster_",sample[i]))

# }

# # Tumor identify 
# # However this not used in tumor classify
# library(Seurat)
# library(pheatmap)
# library(scCNA)
# data(CpRMAP)

# filepath <- "/public/workspace/lily/Lung2Brain/Data/"
# sample <- dir(filepath)[grep("\\.rds$",dir(filepath))] #14

# filepath <- "//public/workspace/lily/Lung2Brain/Data/cluster_data/"
# sample <- c("cluster_BT1299.rds","cluster_BT1300.rds")

# for(i in 1:length(sample)){
# 	# read data
# 	dat <- readRDS(paste0(filepath,sample[i]))

# 	dat@active.ident -> tag
# 	#names(tag) = rownames(info.f)
# 	dat[['RNA']]@data -> dat.f

# 	rs = identificationTumorCells(mat.f=log2(dat.f+1), cluster.tag=tag, chrpq=GRCh37.chrpq)

# 	#get the result information 
# 	save(rs,file=paste0('/public/workspace/lily/Lung2Brain/Data/rs/',gsub("\\.rds$","",sample[i]),"_rs.RData"))
# 	#rs$information$putativeTumor2
# 	dat <- AddMetaData(object = dat, metadata = rs$information$putativeTumor2, col.name = "putativeTumor2")

# 	dat <- AddMetaData(object = dat, metadata = rs$information$putativeTumor, col.name = "putativeTumor")

# 	saveRDS(dat, file = paste0('/public/workspace/lily/Lung2Brain/Data/tumoridentify/',gsub("\\.rds$","",sample[i]),"_identify.rds"))

# }









#####################################################################################################################################################
#####################################################################################################################################################
# integration sample 
#
library(Seurat)
inte.list <- list()
names <- c()
sample <- dir("/public/workspace/lily/Lung2Brain/Data/tumoridentify/")
sample <- sample[c(1,2,6,7,16,17,19)]
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
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte, dims = 1:20)
inte <- FindClusters(inte, resolution = 2)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte,dims = 1:20)
inte$type_group <- "LC"
inte$type_group[which(inte$orig.ident%in%c("A20190305","T-Bsc1","A20190312"))] <- "LCBM"
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/inte7/inte7.RDS")

#===========================================================================================================================================
# cell type annonation
# 2021-4-17

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7.RDS")
DefaultAssay(dat) <- "RNA"
pdf("/public/workspace/lily/Lung2Brain/inte7/inte7_marker_featureplot.pdf")
DimPlot(dat,label=T)
FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,order=T) # T cell
FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','FCGR1A'),label=T,order=T) # Myeloid
FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,order=T) # B cell 
FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,order=T) # Oligodendrocyte
FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,order=T) # Fibroblast/Vascular
FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,order=T) # Endothelial
FeaturePlot(dat,features=c("EGFR","MET","KRAS"),label=T,order=T,cols=c("lightgrey", "red"))
dev.off()

#===========================================================================================================================================
dat$type <- "unknow"
dat$type[which(dat$seurat_clusters%in%c(1,2,3,18,39,36,9,30,19,34,40,43,15,28,35,38))] <- "T_cell"
dat$type[which(dat$seurat_clusters%in%c(4,6,12,13,16,7,26,25,14))] <- "Myeloid"
dat$type[which(dat$seurat_clusters%in%c(11,23))] <- "B_cell"
dat$type[which(dat$seurat_clusters%in%c(29,27))] <- "Oligodendrocyte"
dat$type[which(dat$seurat_clusters%in%c(22,33))] <- "Fibroblast"
dat$type[which(dat$seurat_clusters%in%c(32))] <- "Endothelial"
dat$type[which(dat$seurat_clusters%in%c(0,5,10,20,21,37,24,17,8,31))] <- "maliganant"


#==========================================================================================================================================
dat$maliganant <- "non-tumor"
dat$maliganant[which(dat$type=="maliganant")] <- "tumor"

data <- dat[['RNA']]@data
mean(apply(data[,which(dat$maliganant=="non-tumor")],2,function(x){length(which(x[]>0))}))

mean(apply(data[,which(dat$maliganant=="tumor")],2,function(x){length(which(x[]>0))}))

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")







######################################################################################################################################################
######################################################################################################################################################
# integration 9 sample to analysise TME 
# 2021-4-17 record
# this way integrate 9 sample is OK 
#=====================================================================================================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")
GBM1 <- readRDS("/public/workspace/lily/PS/Final_716/lesion1.RDS")
GBM2 <- readRDS("/public/workspace/lily/PS/Final_716/lesion2.RDS")

integration.anchors <- FindIntegrationAnchors(object.list = c(dat,GBM1,GBM2))
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

#=================================================================================================================================
# add some information 
#
#=================================================================================================================================
inte$type_group[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- "GBM" # sample info 

# cell type info 
inte$type[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- inte$celltype[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))]
inte$type[which(inte$type%in%c("BMDM","MG"))] <- "Myeloid"
inte$type[which(inte$type=="Tumor_cell")] <- "maliganant"

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
















# InferCNV 
#================================================================================================================================================
library(infercnv)
#================================================================================================================================================

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")
sample <- sample <- unique(dat$orig.ident)
# celltype <- unique(dat$type)
# celltype <- celltype[-grep("unknow",celltype)]

for(i in 1:length(sample)){
	tmp <- subset(dat,cells=which(dat$orig.ident==sample[i]&dat$type%in%c("T_cell","maliganant")))
	dir.create(paste0("/public/workspace/lily/Lung2Brain/inte7/infercnv/",sample[i]))
	setwd(paste0("/public/workspace/lily/Lung2Brain/inte7/infercnv/",sample[i]))

	write.table(tmp$type,paste(sample[i],"cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
	count<-as.matrix(tmp@assays$RNA@counts)
	write.table(count,paste(sample[i],"count_exp.txt",sep='_'),sep="\t",quote=F)

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(sample[i],"count_exp.txt",sep='_'),
         annotations_file=paste(sample[i],"cell_info.txt",sep='_'),
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names=c("T_cell"))
         

	infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=10, #big
                             no_plot=F ,
                             denoise=TRUE,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )


	plot_cnv(infercnv_obj, out_dir = "./",output_filename = paste0("infercnv",sample[i]), output_format = "pdf",color_safe_pal=F)


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE, 
                             ,
                             HMM=TRUE)

}





































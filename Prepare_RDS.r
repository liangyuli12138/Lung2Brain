
# code record 
# this code in 2021-12-30
#==========================================================================================================================================
# 1. run cellranger 
# cell ranger was run by hg38,use cellranger v3.0.2
# cell ranger was run by xiapeng and lly
# xiapeng run all LCBM samples and lly run paired Lung sample
#==========================================================================================================================================
# 2. sample filter 
# Lung cancer data [Sangsung] were filter by MT<20%.
# sample for Sangsung was prepare in Prepare_GSE131907.R file 
# so maybe LCBM also use MT <20%
############################################################################################################################################
############################################################################################################################################
#===========================================================================================================================================
# @ this sample do not need to normalized 
# LCBM cell ranger result in xiapeng's file  /public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/
sample <- grep("\\.|E20927",dir("/public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/"),invert=T,value=T)

for(i in 1:length(sample)){
	name <- sample[i]
	datafile = paste0("/public/workspace/xiapeng/Brain_Tumor_sc/0data/LCBM/inhouse/",name,"/3cellranger/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/HG38_Data/Prepare/"
	rdspath = '/public/workspace/lily/Lung2Brain/HG38_Data/RDS/'

	tmp <- Read10X(data.dir = datafile)

    if(name=="WSY"){
        name <- "Pair_BM"
    }
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()

    if(ncol(dat)>10000){
        dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10) 
    }else{
        dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20) 
    }
	
   	dat = NormalizeData(object = dat)
	dat <- FindVariableFeatures(object = dat)
	saveRDS(dat,file=paste0(rdspath,name,'.RDS'))
}






#===========================================================================================================================================
#  and for paired Lung samples
# ~/R21136163_hg38 ,will mv into /public/workspace/lily/LCBM_pair/

library(Seurat)
	name <- "Pair_Lung"
	datafile = paste0("/public/workspace/lily/","R21136163_hg38","/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/HG38_Data/Prepare/"
	rdspath = '/public/workspace/lily/Lung2Brain/HG38_Data/RDS/'

	tmp <- Read10X(data.dir = datafile)
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()

	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
   	dat = NormalizeData(object = dat)
	dat <- FindVariableFeatures(object = dat)
	saveRDS(dat,file=paste0(rdspath,name,'.RDS'))


############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
# change some name 
# mv 
# and this data will trans to .5 
# scp lily@202.195.187.3:/public/workspace/lily/Lung2Brain/HG38_Data/RDS/\{A20190305.RDS,A20190312.RDS,D0927.RDS,E0927.RDS,Pair_BM.RDS,Pair_LUNG.RDS,T_Bsc1.RDS\} ./











# 2022-2-14
# prepare RDS for 
#================================================================================================================================
library(Seurat)
	name <- "PLCBM2"
	datafile = paste0("/public/workspace/lily/LCBM_pair/","R22009109_hg38","/outs/filtered_feature_bc_matrix")
	outpath = "/public/workspace/lily/Lung2Brain/HG38_Data/Prepare/"
	rdspath = '/public/workspace/lily/Lung2Brain/HG38_Data/RDS/'

	tmp <- Read10X(data.dir = datafile)
	dat<- CreateSeuratObject(counts = tmp,  project = name,min.cells = 3, min.features = 200)
	dat[["percent.mt"]] <- PercentageFeatureSet(object = dat, pattern = "^MT-")
	pdf(paste0(outpath,name, "_vlnPlot_prepare.pdf"))
	p<-VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	print(p)
	dev.off()

	dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
   	dat = NormalizeData(object = dat)
	dat <- FindVariableFeatures(object = dat)
	saveRDS(dat,file=paste0(rdspath,name,'.RDS'))

# and this data will trans to .5 
# scp lily@202.195.187.3:/public/workspace/lily/Lung2Brain/HG38_Data/RDS/PLCBM2.RDS ./













































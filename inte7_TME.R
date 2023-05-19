#!/usr/bin/Rscript
#==============================================================================================
# 2020-8-25
# Lung2Brain 
# lly
# inte7 TME analysize 
#
#==============================================================================================


#=================================================================================================================================
# TME analysis 
# 
#=================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
dat$malignant[which(dat$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271') & dat$type=="malignant")] <- "malignant"
#======================= subset cells ======================================
ntumor <- subset(dat,cells=which(dat$malignant=="non-tumor"))

#=================================================================================================================================
library(Seurat)
ntumor <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_ntumor_recluster.RDS")
table(ntumor$type,ntumor$type_group) -> tmp
apply(tmp,2,function(x){x[]/sum(x)}) -> a
barplot(a,col=c('#377EB8','#910241','#984EA3','#F29403','#8CA77B','#B2DF8A','#999999'))

cols <- c('#377EB8','#910241','#984EA3','#E41A1C','#F29403','#8CA77B','#B2DF8A','#999999')

DimPlot(ntumor,group.by="type",cols=c('#377EB8','#910241','#984EA3','#F29403','#8CA77B','#B2DF8A','#999999'))






#==============================================================================================================================
# Myeloid cell 
# BMDM and MG 
#==============================================================================================================================
library(Seurat)
ntumor <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_ntumor_recluster.RDS")
tmp <- subset(ntumor,cells=which(ntumor$type=="Myeloid"))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tmp[['RNA']]@data),c("BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- as.data.frame(mod)
#==================================================================
# add cell type annonation 
#
#==================================================================

apply(mod[,5:6],1,function(x){
	tmp <- length(names(which(x[]<0.05)))
	if(tmp==0|tmp>1){
		c("NA")
	}else{
		names(which(x[]<0.05))
	}
}) -> mod$type

mod$type <- gsub("_marker_pval","",as.vector(mod$type))
mod$type[which(mod$type=="NA")] <- "Myeloid"
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/Myeloid_type.RDS")





#==============================================================================================================================
# T cell 
# CD4/CD8 Treg 
#==============================================================================================================================
library(Seurat)
ntumor <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_ntumor_recluster.RDS")
tdat <- subset(ntumor,cells=which(ntumor$type=="T_cell"))
recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:20)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=1)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:20,check_duplicates = FALSE)

	return(tmp_dat)
}
tdat <- recluster(tdat)
saveRDS(tdat,file="/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.RDS")

#======================================================================================
# plot and classify t cell subtype 
#
#======================================================================================
tdat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.RDS")

#======================================================================================
# filter some cells
# use CD3D or CD3E
#======================================================================================
data <- tdat[['RNA']]@data
tdat.f <- subset(tdat,cells=names(which(data["CD3D",]>0|data["CD3E",]>0)))

#=====================================================================================
# T cell need to re-intergrated 
#
#=====================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")
inte.list <- list() 
samplelist <- unique(dat$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}
#====================================================================================
# integrated t cell 
# 

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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")


#======================================================================================
# plot 
#======================================================================================
tdat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")
DefaultAssay(tdat) <- "RNA"
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Tcell_featureplot.pdf",useDingbats=F)
DimPlot(tdat,label=T)
DimPlot(tdat,group.by="type_group")
FeaturePlot(tdat,features=c("CD4","IL7R","CD3D"),label=T,order=T) # CD4 T
FeaturePlot(tdat,features=c("CD8A","CD8B","CD3D"),label=T,order=T) # CD8 T 
FeaturePlot(tdat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),label=F,order=T) # Treg 
FeaturePlot(tdat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),label=F,order=T) # exhausted
# FeaturePlot(tdat,features=c("TCF7","SELL","LEF1","CCR7"),label=T,order=T) # naive T cell  
# FeaturePlot(tdat,features=c("IL2","GZMA","GNLY"),label=T,order=T) # cytotoxic
dev.off()

  






#================================================================================================================================
# use GSE131907 to analysis 
#================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE131907/cell_ann.RDS")
dat.f <- subset(dat,cells=which(dat$Cell_type.refined=="T/NK cells"))
dat.ff <- subset(dat.f,cells=which(dat.f$Sample_Origin%in%c("nLung","tL/B","tLung")))
data <- dat.ff[['RNA']]@data
dat.fff <- subset(dat.ff,cells=names(which(data["CD3D",]>0|data["CD3E",]>0)))
saveRDS(dat.fff,file="/public/workspace/lily/metastasis/data/verify/GSE131907/T_cell.RDS")



# 2020-10-3 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE131907/T_cell.RDS")
tmp <- subset(dat,cells=which(dat$Sample_Origin%in%c("tL/B","tLung")))
tmp <- recluster(tmp)

res <- apply(table(dat$Cell_subtype,dat$Sample_Origin),2,function(x){x[]/sum(x)})
tmp <- table(dat$Cell_subtype,dat$Sample_Origin)

pdf("./GSE131907_Tcell_percentage.pdf")
barplot(res["Cytotoxic CD8+ T",c("tL/B","tLung")],col=c("red","blue"),main="Cytotoxic CD8+ T")
barplot(res["Exhausted CD8+ T",c("tL/B","tLung")],col=c("red","blue"),main="Exhausted CD8+ T")
dev.off()


#===============================
# cell annoation
#===============================
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE131907_Tcell_marker.pdf")
DimPlot(tmp,group.by="Sample_Origin")
FeaturePlot(tmp,features=c("IL7R","CD4","CD8A","CD8B"))
FeaturePlot(tmp,features=c("LAG3", "TIGIT", "PDCD1", "HAVCR2"))
dev.off()




#==============================================================================================================================
# cellphone DB for each sample 
# /public/workspace/lily/Lung2Brain/TME/CCC_data/
#==============================================================================================================================

library(Seurat)
cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}


dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS") 
# have brain metastasis,lung adenocarcinoma and glioblastoma data
samplelist <- as.vector(unique(dat$orig.ident))
for(i in 1:length(samplelist)){
	tmp <-subset(dat,cells=which(dat$orig.ident==samplelist[i]))
	dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",samplelist[i]))
	outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",samplelist[i],"/")
	cellphoneDB_input(tmp,"type",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))
}


#================================================================================================
# run in bash 
# 
#================================================================================================


for i in `ls /public/workspace/lily/Lung2Brain/TME/CCC_data`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/CCC_data/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/CCC_data/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/CCC_data/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/CCC_data/${i}/expr.txt --threads 10 \
	--output-path=${RESPATH} --counts-data gene_name

done







#=============================================================================================================================================================================================
# check cellphoneDB result 
#
#=============================================================================================================================================================================================
# T cell and tumor cell interaction analysis 
#================================================================
LUAD <- c("scrBT1432m","scrBT1431m","BT1296","BT1297")
LCBM <- c("A20190305","A20190312","T_Bsc1")
GBM <- c("lesion1","lesion2")

#================================================================
# lly have test if the interaction pair is unique 
# 2020-8-25
#===============================
for(i in 1:length(LUAD)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",LUAD[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",LUAD[i],"/res/pvalues.txt"),sep="\t",header=T)

#======== change names
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair

#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(pvalue))]
	means.f <- means[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",LUAD[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",LUAD[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(LUAD[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(BT1296[[1]]), # pvalue data 
 					rownames(BT1297[[1]]),
 					rownames(scrBT1431m[[1]]),
 					rownames(scrBT1432m[[1]])
 					))

copvalue <-cbind(BT1297[[1]][copair,],BT1296[[1]][copair,],scrBT1431m[[1]][copair,],scrBT1432m[[1]][copair,])
copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]

#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
ann <- data.frame(row.names=colnames(copvalue.f),interact_type=gsub("_BT.*|_scr.*","",colnames(copvalue.f)),sample=sapply(strsplit(colnames(copvalue.f),"_"),function(x){x[length(x)]}))
ann<- ann[order(ann$interact_type),]

library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LUAD_copvalue_heatmap.pdf",height=15)
pheatmap(copvalue.f[,rownames(ann)],cluster_col=F,annotation_col = ann,color = colorRampPalette(c("steelblue",'white','red'))(500),border=F)
dev.off()

#==================================================================
# save Data
# 
#==================================================================
saveRDS(copvalue.f,file="/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LUAD_copvalue.f.RDS")



















#==========================================================================================================================================
# LCBM 
# change T-Bsc1 to T_Bsc1
#==========================================================================================================================================
LCBM <- c("A20190305","A20190312","T_Bsc1")

for(i in 1:length(LCBM)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",LCBM[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",LCBM[i],"/res/pvalues.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair

#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(pvalue))]
	means.f <- means[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",LCBM[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",LCBM[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(LCBM[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(A20190305[[1]]), # pvalue data 
 					rownames(A20190312[[1]]),
 					rownames(T_Bsc1[[1]])
 					))

copvalue <-cbind(A20190305[[1]][copair,],A20190312[[1]][copair,],T_Bsc1[[1]][copair,])
copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]

#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
ann <- data.frame(row.names=colnames(copvalue.f),interact_type=gsub("_A2019.*|_T_Bsc.*","",colnames(copvalue.f)),sample=sapply(strsplit(colnames(copvalue.f),"_"),function(x){x[length(x)]}))
ann<- ann[order(ann$interact_type),]

library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LUAD_copvalue_heatmap.pdf",height=15)
pheatmap(copvalue.f[,rownames(ann)],cluster_col=F,annotation_col = ann,color = colorRampPalette(c("steelblue",'white','red'))(500),border=F)
dev.off()


#==================================================================
# save Data
# 
#==================================================================
saveRDS(copvalue.f,file="/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LCBM_copvalue.f.RDS")















#==========================================================================================================================================
# GBM 
# lesion1 do not have T_cell, so maybe need to change a way to show 
#==========================================================================================================================================
GBM <- c("lesion1","lesion2")

for(i in 1:length(GBM)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",GBM[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",GBM[i],"/res/pvalues.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair

#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(pvalue)),drop=F]
	means.f <- means[,grep("malignant\\.T_cell|malignant\\.Myeloid|Myeloid\\.T_cell",colnames(means)),drop=F]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",GBM[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",GBM[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(GBM[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(lesion1[[1]]), # pvalue data 
 					rownames(lesion2[[1]])
 					))

copvalue <-cbind(lesion1[[1]][copair,],lesion2[[1]][copair,])
copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
colnames(copvalue.f)[1] <- "malignant.Myeloid_lesion1"

#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
ann <- data.frame(row.names=colnames(copvalue.f),interact_type=gsub("_lesion.*","",colnames(copvalue.f)),sample=sapply(strsplit(colnames(copvalue.f),"_"),function(x){x[length(x)]}))
ann<- ann[order(ann$interact_type),]

library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LUAD_copvalue_heatmap.pdf",height=15)
pheatmap(copvalue.f[,rownames(ann)],cluster_col=F,annotation_col = ann,color = colorRampPalette(c("steelblue",'white','red'))(500),border=F)
dev.off()


#==================================================================
# save Data
# 
#==================================================================
saveRDS(copvalue.f,file="/public/workspace/lily/Lung2Brain/TME/cellphoneDB/GBM_copvalue.f.RDS")


















#=====================================================================================================================================================================================================
# malignant to Meyloid cell 
# 2020-8-25
# try to another way to show 
#=====================================================================================================================================================================================================
samplelist <- c("scrBT1432m","scrBT1431m","BT1296","BT1297","A20190305","A20190312","T_Bsc1","lesion1","lesion2")
for(i in 1:length(samplelist)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",samplelist[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CCC_data/",samplelist[i],"/res/pvalues.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair

#======== grep Tumor and T cell 
	if( length(grep("malignant\\.Myeloid",colnames(pvalue)))>0 ){
		pvalue.f <- pvalue[,grep("malignant\\.Myeloid",colnames(pvalue)),drop=F]
		means.f <- means[,grep("malignant\\.Myeloid",colnames(means)),drop=F]

	#======== process data 
		pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
		pvalue.f[pvalue.f>0.05] <- 1
		pvalue.ff <- -log2(pvalue.f)

	#======== change colnames for next merge data 
		colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",samplelist[i])
		colnames(means.f) <- paste0(colnames(means.f),"_",samplelist[i])

	#=============================== 
	# just use pvalue or both 
	#
	#================================
		assign(samplelist[i],list(pvalue.ff,means.f))

	}
	
}




copair <- Reduce(intersect,list( rownames(scrBT1431m[[1]]),rownames(scrBT1432m[[1]]),rownames(A20190305[[1]]),rownames(A20190312[[1]]),rownames(T_Bsc1[[1]]),
					rownames(BT1296[[1]]),rownames(BT1297[[1]]),rownames(lesion1[[1]]),rownames(lesion2[[1]])
 					)
				)

copvalue <-cbind(A20190305[[1]][copair,,drop=F],A20190312[[1]][copair,,drop=F],T_Bsc1[[1]][copair,,drop=F],
				scrBT1431m[[1]][copair,,drop=F],scrBT1432m[[1]][copair,,drop=F],BT1296[[1]][copair,,drop=F],BT1297[[1]][copair,,drop=F],
				lesion1[[1]][copair,,drop=F],lesion2[[1]][copair,,drop=F]
				)
copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]


#===============================================================================================================================
# 
ann <- data.frame(row.names=colnames(copvalue.f),interact_type= substr(colnames(copvalue.f),1,17),sample= substr(colnames(copvalue.f),19,stop=60))

#============================ add some type information
ann$type <- "LUAD"
ann$type[which(ann$sample%in%c("A20190305","A20190312","T_Bsc1"))] <- "LCBM"
ann$type[which(ann$sample%in%c("lesion1","lesion2"))] <- "GBM"
ann <- ann[order(ann$type),] # order by cancer type 

pheatmap(copvalue.f[,rownames(ann)],cluster_col=F,annotation_col = ann,color = colorRampPalette(c("steelblue",'white','red'))(500),border=F)







#====================================================================================================================================================================================
# another way to show 
#
#====================================================================================================================================================================================
# LCBM
#
#===============================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LCBM_copvalue.f.RDS")
res <- dat[,grep("malignant.Myeloid",colnames(dat))]
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_CCC.pdf",height=15)
pheatmap(res[which(rowSums(res)>0),])
dev.off()



dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/LUAD_copvalue.f.RDS")
res <- dat[,grep("malignant.Myeloid",colnames(dat))]
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LUAD_CCC.pdf",height=15)
pheatmap(res[which(rowSums(res)>0),])
dev.off()



dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/cellphoneDB/GBM_copvalue.f.RDS")
res <- dat[,grep("malignant.Myeloid",colnames(dat))]
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GBM_CCC.pdf",height=15)
pheatmap(res[which(rowSums(res)>0),])
dev.off()





















#==========================================================================================================
# exhausted T cell 
#
#==========================================================================================================
library(Seurat)
tdat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")
DefaultAssay(tdat) <- "RNA"
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Tcell_exhausted_featureplot.pdf",useDingbats=F)
DimPlot(tdat,label=T)
DimPlot(tdat,group.by="type_group")
# FeaturePlot(tdat,features=c("CD4","IL7R","CD3D"),label=T,order=T) # CD4 T
# FeaturePlot(tdat,features=c("CD8A","CD8B","CD3D"),label=T,order=T) # CD8 T 
#FeaturePlot(tdat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),label=F,order=T) # Treg 
FeaturePlot(tdat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),label=F,order=T) # exhausted
# FeaturePlot(tdat,features=c("TCF7","SELL","LEF1","CCR7"),label=T,order=T) # naive T cell  
# FeaturePlot(tdat,features=c("IL2","GZMA","GNLY"),label=T,order=T) # cytotoxic
dev.off()








#==========================================================================================================
# 2020-10-14 
# Cellphone DB analysis 
#==========================================================================================================
bytlib load R-3.6.0
R

# cellphone DB result have done 
# BMDM MG different 

#===========================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Mye_type <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/Myeloid_type.RDS")

# purity Myeloid cells
dat$type.refine <- dat$type
dat$type.refine[which(rownames(dat@meta.data)%in%rownames(Mye_type))] <- Mye_type$type

# purity T cells
tcell  <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")
tcell$type.refine <- paste0(tcell$type,"_refine")
dat$type.refine[which(rownames(dat@meta.data)%in%rownames(tcell@meta.data))] <- tcell$type.refine
#=============================================================================================
# just calculate BMDM Tcell MG and tumor interaction
#
#=============================================================================================

tmp <- subset(dat,cells=which(dat$type.refine%in%c("malignant","MG","BMDM","T_cell_refine")))
saveRDS(tmp,file="/public/workspace/lily/Lung2Brain/inte7/Data/TME9_T_MDM_MG.RDS")





#============================================================================================
# run CellphoneDB
#

library(Seurat)
cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}


dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/TME9_T_MDM_MG.RDS") 
# have brain metastasis,lung adenocarcinoma and glioblastoma data
samplelist <- as.vector(unique(dat$orig.ident))
for(i in 1:length(samplelist)){
	tmp <-subset(dat,cells=which(dat$orig.ident==samplelist[i]))
	dir.create(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",samplelist[i]))
	outpath <- paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",samplelist[i],"/")
	cellphoneDB_input(tmp,"type.refine",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))
}


#================================================================================================
# run in bash 
# 
#================================================================================================


for i in `ls /public/workspace/lily/Lung2Brain/inte7/CellphoneDB/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/inte7/CellphoneDB/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/inte7/CellphoneDB/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/inte7/CellphoneDB/${i}/expr.txt --threads 10 \
	--output-path=${RESPATH} --counts-data gene_name

done




#========================================================================================
# check result 
# first check Tumor and BMDM MG result 
#========================================================================================
# BMDM with tumor
#========================================================================================
# BM 
BM_file <- c("A20190305","A20190312","T_Bsc1")
LC_file <- c("BT1296","BT1297","scrBT1431m","scrBT1432m")
GBM_file <- c("lesion2","lesion1")


#========================
#===============================
for(i in 1:length(BM_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",BM_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",BM_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(BM_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(A20190305[[1]]), # pvalue data 
 					rownames(A20190312[[1]]),
 					rownames(T_Bsc1[[1]])
 					))

copvalue <-cbind(A20190305[[1]][copair,],A20190312[[1]][copair,],T_Bsc1[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.MG",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.BM
#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_T2MG_copvalue.pdf",height=20)
pheatmap(copvalue.ff,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()










#==================================================================================================
# LC file 
#
#==================================================================================================

BM_file <- c("A20190305","A20190312","T_Bsc1")
LC_file <- c("BT1296","BT1297","scrBT1431m","scrBT1432m")
GBM_file <- c("lesion2","lesion1")


#========================
#===============================
for(i in 1:length(LC_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",LC_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",LC_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(LC_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(BT1296[[1]]), # pvalue data 
 					rownames(BT1297[[1]]),
 					rownames(scrBT1431m[[1]]),
 					rownames(scrBT1432m[[1]])
 					))

copvalue <-cbind(BT1296[[1]][copair,],BT1297[[1]][copair,],scrBT1431m[[1]][copair,],scrBT1432m[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.BMDM",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.LC
#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/LC_BMDM2T_copvalue.pdf",height=20)
pheatmap(copvalue.ff,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()




merge(copvalue.ff.BM,copvalue.ff.LC,by="row.names") -> res
rownames(res) <- res$Row.names
res$Row.names <- NULL
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/LC_BM_BMDM2T_copvalue.pdf",height=20)
pheatmap(res,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_T2BMDM_unique.pdf",height=20)
pheatmap(copvalue.ff.BM[-which(rownames(copvalue.ff.BM)%in%rownames(res)),],cluster_col=F,color = colorRampPalette(c('white','red'))(500))
dev.off()


#=========================================================================================
# GBM file 
# compare MG with BM 
#=========================================================================================

GBM_file <- c("lesion2","lesion1")

#========================
#===============================
for(i in 1:length(GBM_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",GBM_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",GBM_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",GBM_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",GBM_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(GBM_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(lesion1[[1]]), # pvalue data 
 					rownames(lesion1[[1]])
 					))

copvalue <-cbind(lesion1[[1]][copair,],lesion2[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.MG",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.GBM
#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/GBM_T2MG_copvalue.pdf",height=20)
pheatmap(copvalue.ff,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()



merge(copvalue.ff.BM,copvalue.ff.GBM,by="row.names") -> res
rownames(res) <- res$Row.names
res$Row.names <- NULL
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/GBM_BM_T2MG_copvalue.pdf",height=20)
pheatmap(res,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_T2MG_unique.pdf",height=20)
pheatmap(copvalue.ff.BM[-which(rownames(copvalue.ff.BM)%in%rownames(res)),],cluster_col=F,color = colorRampPalette(c('white','red'))(500))
dev.off()






#==========================================================================================
# T cell analysis 
#
#==========================================================================================
# BM samples 
# cellchat plot 
#==========================================================================================
# BM 
BM_file <- c("A20190305","A20190312","T_Bsc1")

#===============================
#===============================
for(i in 1:length(BM_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG|T_cell_refine",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG|T_cell_refine",colnames(means))]


#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",BM_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",BM_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(BM_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(A20190305[[1]]), # pvalue data 
 					rownames(A20190312[[1]]),
 					rownames(T_Bsc1[[1]])
 					))

copvalue <-cbind(A20190305[[1]][copair,],A20190312[[1]][copair,],T_Bsc1[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.MG",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.BM

# tmp_res <- matrix(nrow=nrow(copvalue),ncol=16)
# for(i in 1:16){
# 	tmp<- mean(copvalue[,i],copvalue[,i+16],copvalue[,i+32])
# }


tmp_res <- apply(copvalue,1,function(x){
	tmp <- c()
	for(i in 1:16){
		tmp <- c(tmp,mean(x[i],x[i+16],x[i+32]))
	}
	tmp
})

rownames(tmp_res) <- gsub("_A20190305$","",colnames(copvalue)[1:16])
res.f <- t(tmp_res)

#===========================================================

stat1 <- c()
for(j in 1:ncol(res.f)){
	tmp_count1 <- 0
	for(i in 1:nrow(res.f)){
		if(res.f[i,j]>0){
			tmp_count1 <- tmp_count1+1
		}	
	}
	stat1 <- c(stat1,tmp_count1)
	names(stat1)[j] <- colnames(res.f)[j]
}

res.m <- matrix(ncol=4)
for(i in 1:4){
	tmp <- stat1[(i-1)*4+(1:4)]
	res.m <- rbind(res.m,tmp)
}
res.m <- res.m[-1,]
colnames(res.m) <- gsub("^BMDM\\.","",colnames(res.m))
rownames(res.m) <- colnames(res.m)

rownames(res.m)[3] <- "T cell"
colnames(res.m)[3] <- "T cell"


.libPaths() -> tmp
.libPaths(c(tmp,"/public/workspace/zhumy/R/x86_64-pc-linux-gnu-library/3.6.0"))
library(CellChat)

res.mf <- apply(res.m,2,function(x){scale(x,center=F)})
rownames(res.mf) <- colnames(res.mf)
cairo_pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_CellChat.pdf")
# Myeloid Oligodendrocyte  T_cell Tumor_cell
cols <- c('#F29403','#A6CEE3','#4DAF4A',"#E41A1C")
netVisual_circle(res.m,weight.scale=T,color.use=cols)
dev.off()



#====================================================================================================
# LC cellchat 
#
#====================================================================================================
LC_file <- c("BT1296","BT1297","scrBT1431m","scrBT1432m")

#===============================
#===============================
for(i in 1:length(LC_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue <- pvalue[,12:ncol(pvalue)]
	means <- means[,12:ncol(means)]
	if(length(grep("MG",colnames(pvalue)))>0){
		pvalue.f <- pvalue[,-grep("MG",colnames(pvalue))]
		means.f <- means[,-grep("MG",colnames(means))]
	}else{
		pvalue.f <- pvalue
		means.f <- means
	}
	
	


#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",LC_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",LC_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(LC_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(BT1296[[1]]), # pvalue data 
 					rownames(BT1297[[1]]),
 					rownames(scrBT1431m[[1]]),
 					rownames(scrBT1432m[[1]])
 					))

copvalue <-cbind(BT1296[[1]][copair,],BT1297[[1]][copair,],scrBT1431m[[1]][copair,],scrBT1432m[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.BMDM",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.LC

tmp_res <- apply(copvalue,1,function(x){
	tmp <- c()
	for(i in 1:9){
		tmp <- c(tmp,mean(x[i],x[i+9],x[i+18],x[i+27]))
	}
	tmp
})

rownames(tmp_res) <- gsub("_BT1296$","",colnames(copvalue)[1:9])
res.f <- t(tmp_res)

#===========================================================

stat1 <- c()
for(j in 1:ncol(res.f)){
	tmp_count1 <- 0
	for(i in 1:nrow(res.f)){
		if(res.f[i,j]>0){
			tmp_count1 <- tmp_count1+1
		}	
	}
	stat1 <- c(stat1,tmp_count1)
	names(stat1)[j] <- colnames(res.f)[j]
}

res.m <- matrix(ncol=3)
for(i in 1:3){
	tmp <- stat1[(i-1)*3+(1:3)]
	res.m <- rbind(res.m,tmp)
}
res.m <- res.m[-1,]
colnames(res.m) <- gsub("^BMDM\\.","",colnames(res.m))
rownames(res.m) <- colnames(res.m)

rownames(res.m)[2] <- "T cell"
colnames(res.m)[2] <- "T cell"


.libPaths() -> tmp
.libPaths(c(tmp,"/public/workspace/zhumy/R/x86_64-pc-linux-gnu-library/3.6.0"))
library(CellChat)

res.mf <- apply(res.m,2,function(x){scale(x,center=F)})
rownames(res.mf) <- colnames(res.mf)
cairo_pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/LC_CellChat.pdf")
# Myeloid Oligodendrocyte  T_cell Tumor_cell
cols <- c('#F29403','#4DAF4A',"#E41A1C")
netVisual_circle(res.m,weight.scale=T,color.use=cols)
dev.off()



#==========================================================================================
# cellphone DB result 
#
#==========================================================================================

# BM 
BM_file <- c("A20190305","A20190312","T_Bsc1")

#===============================
#===============================
for(i in 1:length(BM_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",BM_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG|T_cell_refine",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG|T_cell_refine",colnames(means))]


#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",BM_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",BM_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(BM_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(A20190305[[1]]), # pvalue data 
 					rownames(A20190312[[1]]),
 					rownames(T_Bsc1[[1]])
 					))

copvalue <-cbind(A20190305[[1]][copair,],A20190312[[1]][copair,],T_Bsc1[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.BMDM",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.BM


library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_T2tcell_copvalue.pdf",height=20)
pheatmap(copvalue.ff,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()



#===========================================================================================
# LC 
#===========================================================================================

BM_file <- c("A20190305","A20190312","T_Bsc1")
LC_file <- c("BT1296","BT1297","scrBT1431m","scrBT1432m")
GBM_file <- c("lesion2","lesion1")


#========================
#===============================
for(i in 1:length(LC_file)){
	pvalue <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/pvalues.txt"),sep="\t",header=T)
	means <- read.table(paste0("/public/workspace/lily/Lung2Brain/inte7/CellphoneDB/",LC_file[i],"/res/means.txt"),sep="\t",header=T)

#======== change names
	if(length(grep("CALCA_CALCR",pvalue$interacting_pair))>1){
		pvalue <- pvalue[-grep("CALCA_CALCR",pvalue$interacting_pair)[2],]
		means <- means[-grep("CALCA_CALCR",means$interacting_pair)[2],]
	}
	rownames(pvalue) <- pvalue$interacting_pair
	rownames(means) <- means$interacting_pair


#======== grep Tumor and T cell 
	pvalue.f <- pvalue[,grep("malignant|BMDM|MG",colnames(pvalue))]
	means.f <- means[,grep("malignant|BMDM|MG",colnames(means))]

#======== process data 
	pvalue.f[pvalue.f==0] <- min(pvalue.f[pvalue.f>0])
	pvalue.f[pvalue.f>0.05] <- 1
	pvalue.ff <- -log2(pvalue.f)

#======== change colnames for next merge data 
	colnames(pvalue.ff) <- paste0(colnames(pvalue.ff),"_",LC_file[i])
	colnames(means.f) <- paste0(colnames(means.f),"_",LC_file[i])

#=============================== 
# just use pvalue or both 
#
#================================
	assign(LC_file[i],list(pvalue.ff,means.f))

}

#================================================================== 
# intersect of one cancer type 
copair <- Reduce(intersect,list(
					rownames(BT1296[[1]]), # pvalue data 
 					rownames(BT1297[[1]]),
 					rownames(scrBT1431m[[1]]),
 					rownames(scrBT1432m[[1]])
 					))

copvalue <-cbind(BT1296[[1]][copair,],BT1297[[1]][copair,],scrBT1431m[[1]][copair,],scrBT1432m[[1]][copair,])
#copvalue.f <- copvalue[-which(rowSums(copvalue)==0),]
copvalue.f <- copvalue[,grep("malignant\\.BMDM",colnames(copvalue))]
copvalue.ff <- copvalue.f[-which(rowSums(copvalue.f)==0),]
copvalue.ff -> copvalue.ff.LC
#==================================================================
# prepare for pheatmap 
# add some sample info 
#==================================================================
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/LC_T2tcell_copvalue.pdf",height=20)
pheatmap(copvalue.ff.LC,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()




merge(copvalue.ff.BM,copvalue.ff.LC,by="row.names") -> res
rownames(res) <- res$Row.names
res$Row.names <- NULL
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/LC_BM_T2tcell_copvalue.pdf",height=20)
pheatmap(res,cluster_col=F,color = colorRampPalette(c('white','red'))(500),border=F)
dev.off()


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCC/BM_T2tcell_unique.pdf",height=20)
pheatmap(copvalue.ff.BM[-which(rownames(copvalue.ff.BM)%in%rownames(res)),],cluster_col=F,color = colorRampPalette(c('white','red'))(500))
dev.off()





#==================================================================================
# VN plot to show 
#
#==================================================================================
library(VennDiagram)

venn.diagram(x=list(LUAD=rownames(copvalue.ff.LC),BM=rownames(copvalue.ff.BM)),fill=c("#4DAF4A","#377EB8"),
	filename="LC_BM_T2BMDM.tif",imagetype="tiff")







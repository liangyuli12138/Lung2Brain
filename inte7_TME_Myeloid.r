
#===============================================================================
# Lung2Brain Myeloid analysis 
# 
#===============================================================================



#################################################################################
# should use absolute and immuncell 
# 2021-1-8
#================================================================================
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/LUAD_absolute.txt",sep="\t",header=T)
tmp.dat$SampleID <- substr(gsub("-",".",tmp.dat$Sample.ID),1,15)
tmp.res <- aggregate(ABSOLUTE~SampleID,data=tmp.dat,FUN=median)
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")

rownames(tmp.res) <- tmp.res$SampleID
res.f <- merge(tmp.res,luad_mod,by="row.names")

# ggplot plot result 
#################################################################################
# 2021-1-9
library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_purity_BMS_cor.pdf",useDingbats=F)
p<-ggplot(res.f,aes(x=BMS_test_norm,y=ABSOLUTE)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Purity") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()



#################################################################################
# immune cell proportion analysis
# 2021-1-9
#################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
# plot 
tmp <- table(dat$type_group,dat$type)
tmp.f <- tmp[,-grep("malignant|unknow",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
##################################################################################
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,level=c("GBM","LCBM","LC"))
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Cell_type_inte.pdf",useDingbats=F)
cols <- c('#377EB8','#910241','#984EA3','#F29403','#8CA77B','#B2DF8A')
c("#f99104","#00b7c9")
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()






























# 1. TME analysis should add GBM data 
# 2. 
#===============================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- subset(dat,cells=which(dat$type=="Myeloid"))

inte.list <- list() 
samplelist <- unique(Myeloid$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(Myeloid,cells=which(Myeloid$orig.ident==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
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
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
#==================================================
# RNA assay re-cluster 
# Better do not use re-cluster's clustering information 
#==================================================
# DefaultAssay(Myeloid) <- "integrated"
# Myeloid <- FindVariableFeatures(object = Myeloid)
# # scaling
# all.genes <- rownames(x = Myeloid)
# Myeloid <- ScaleData(object = Myeloid, features = all.genes)
# # PCA
# Myeloid <- RunPCA(object = Myeloid, features = VariableFeatures(object = Myeloid))
# # clustering
# Myeloid <- FindNeighbors(object = Myeloid,dims=1:10)
# # select proper resolution
# Myeloid <- FindClusters(object = Myeloid)
# # T-SNE
# Myeloid <- RunTSNE(object = Myeloid,dims=1:10,check_duplicates = FALSE)
# Myeloid <- RunUMAP(Myeloid,dims=1:10)

#saveRDS(Myeloid,file="/public/workspace/lily/Lung2Brain/TME/MMyeloid/inte_Myeloid.RDS")


#================================================================================
# calculate BMDM and MG difference 
#================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
mod <- mod.analyze2(as.matrix(dat[["RNA"]]@data),c("MG_marker","BMDM_marker"),"/public/workspace/lily/MOD_file/",permN=1000)

#================================================================================
# add some info to reknow BMDM and MG 
#================================================================================
mod <- data.frame(mod)
mod$type.refine <- "unknow"
mod$type.refine[which(mod$BMDM_marker_norm>mod$MG_marker_norm)] <- "BMDM"
mod$type.refine[which(mod$BMDM_marker_norm<mod$MG_marker_norm)] <- "MG"
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_MG_mod.RDS")
dat$type.refine <- mod$type.refine
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")



#================================================================================
# plot BMDM result 
#================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Myeloid_UMAP.pdf",useDingbats=F)
DimPlot(dat)
DimPlot(dat,group.by="type_group",cols=c("#e63946","#2a9d8f","#f4a261"))
DimPlot(dat,group.by="type.refine")
dev.off()



#================================================================================
# plot alluvial plot 
# ggplot2  alluvial diagram
#================================================================================
library(reshape)
library(ggplot2)
library(ggalluvial)
tmp <- table(dat$type_group,dat$type.refine)
tmp.res <- apply(tmp,1,function(x){x/sum(x)})
# barplot(tmp.res)
# by samples is the same result 
# table(dat$orig.ident,dat$type.refine) -> tmp
# tmp.res <- apply(tmp,1,function(x){x/sum(x)})
# colnames(tmp.res)[5] <- "rd001"
# colnames(tmp.res)[6] <- "rd002"

tmp.dat <- melt(tmp.res,id="col.names")
colnames(tmp.dat) <- c("Type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,level=c("GBM","LCBM","LC"))
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_MG_bar.pdf")
cols <- c("#f99104","#00b7c9")
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Type,stratum = Type, alluvium = Type)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()






#===================================================================================================
# TCGA GBM data and GSE14108 to verify BMDM and MG result 
#===================================================================================================
dat1 <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108-GPL570_series_matrix.txt",comment.char="!",sep="\t",header=T)
dat2 <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108-GPL96_series_matrix.txt",comment.char="!",sep="\t",header=T)
# first run in shell to get some info 
gpl96 <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/GPL96_info.txt",sep="\t",header=T)
gpl570 <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/GPL570_info.txt",sep="\t",header=T)

#===================================================================================================
# some prepare 
#===================================================================================================
# GPL570
dat.570 <- merge(dat1,gpl570,by.x="ID_REF",by.y="ID")
dat.570.f <- dat.570[-grep("///",dat.570$Gene.Symbol),]
dat.570.f$ID_REF <- NULL
res.570 <- aggregate(.~Gene.Symbol,data=dat.570.f,FUN=median)
res.570 <- res.570[-1,]
saveRDS(res.570,file="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_GPL570.RDS")

# GPL96
dat.96 <- merge(dat2,gpl96,by.x="ID_REF",by.y="ID")
dat.96.f <- dat.96[-grep("///",dat.96$Gene.Symbol),]
dat.96.f$ID_REF <- NULL
res.96 <- aggregate(.~Gene.Symbol,data=dat.96.f,FUN=median)
res.96 <- res.96[-1,]
saveRDS(res.96,file="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_GPL96.RDS")

# merge two platform data
dat1 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_GPL570.RDS")
dat2 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_GPL96.RDS")
dat.res <- merge(dat1,dat2,by="Gene.Symbol")

rownames(dat.res) <- dat.res$Gene.Symbol
dat.res$Gene.Symbol <- NULL
saveRDS(dat.res,file="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")


# TCGA GBM data 
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/genomicMatrix",sep="\t",header=T)

#=======================================================================================
# calculate BMDM and MG score 
#=======================================================================================
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
dat.tcga <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/genomicMatrix",sep="\t",header=T)
rownames(dat.tcga) <- dat.tcga$sample
dat.tcga$sample <- NULL
dat.res <- merge(dat.BM,dat.tcga,by="row.names")
rownames(dat.res) <- dat.res$Row.names
dat.res$Row.names <- NULL
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.res),c("MG_marker","BMDM_marker"),"/public/workspace/lily/MOD_file/",permN=1000)

#======================================================================================
mod <- data.frame(mod)
mod$type <- "LCBM"
mod$type[grep("TCGA",rownames(mod))] <- "GBM"

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/GBM_LCBM_BMDM_mod.RDS")
#=======================================================================================
library(ggplot2)
library(reshape2)
mod <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/GBM_LCBM_BMDM_mod.RDS")
mod.1 <- mod[which(mod$type=="LCBM"),]
mod.2 <- mod[which(mod$type=="GBM"),]

mod.1$sample <- rownames(mod.1)
dat.1 <- melt(mod.1[,c("MG_marker_norm","BMDM_marker_norm","sample")],id="sample")

mod.2$sample <- rownames(mod.2)
dat.2 <- melt(mod.2[,c("MG_marker_norm","BMDM_marker_norm","sample")],id="sample")



pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/GBM_LCBM_BMDM_new.pdf",useDingbats=F)
ggplot(dat.1, aes(x=value, fill=variable)) + geom_density()+theme_bw()+xlab("LCBM")
ggplot(dat.2, aes(x=value, fill=variable)) + geom_density()+theme_bw()+xlab("GBM")
dev.off()







pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/GBM_LCBM_BMDM.pdf",useDingbats=F)
ggplot(mod, aes(x=MG_marker_norm, y=BMDM_marker_norm, group=type)) + geom_point()

# density plot 
ggplot(mod, aes(x=MG_marker_norm, fill=type)) + geom_density()+theme_bw()

ggplot(mod, aes(x=BMDM_marker_norm, fill=type)) + geom_density()+theme_bw()
dev.off()



#======================================================================================
# Myeloid to different M1 and M2 myeloid 
# the gene marker is come from xiaoxuej 
# just use add module ?
#======================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
subset.dat <- subset(dat,cells=which(dat$type.refine=="BMDM"))
M1_marker <- c("FCGR1A","FCGR1B","CD86","CXCL10","GBP1","CXCL9","IL1B","IL6","CXCL11","TNF","IL23A")
M2_marker <- c("FCER2","IL27RA","CLEC4A","CCL22","CCL18","CCL17","F13A1","IL10","MRC1","CD163")

# use ssgsea to calculate 
#======================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(M1_marker,'M1_marker',out='/public/workspace/lily/MOD_file/M1_marker.mod')
mod.generate(M2_marker,'M2_marker',out='/public/workspace/lily/MOD_file/M2_marker.mod')

mod <- mod.analyze2(as.matrix(subset.dat[['RNA']]@data),c("M1_marker","M2_marker"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_M1_M2.RDS")


# subset.dat <- AddModuleScore(object = subset.dat,assay="RNA",features = M1_marker,name = 'M1_feature')
# subset.dat <- AddModuleScore(object = subset.dat,assay="RNA",features = M2_marker,name = 'M2_feature')

#======================================================================================
# check result 
# calculated result show that M1 and M2 in GBM is not ok 
#======================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
subset.dat <- subset(dat,cells=which(dat$type.refine=="BMDM"))

mod <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_M1_M2.RDS")
mod <- data.frame(mod)
mod$M1_M2_type <- "unknow"
mod$M1_M2_type[mod$M2_marker_norm > mod$M1_marker_norm] <- "M2"
mod$M1_M2_type[mod$M2_marker_norm < mod$M1_marker_norm] <- "M1"
subset.dat$M1_M2_type <- mod$M1_M2_type




















#===============================================================================
# use tumor cell and BMDM/MG to run CellphoneDB 
#===============================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type%in%c("malignant","Myeloid")))
sub.dat$type.refine <- sub.dat$type 
sub.dat$type.refine[which(colnames(sub.dat)%in%colnames(Myeloid))] <- Myeloid$type.refine
sub.dat$sample <- sub.dat$orig.ident
sub.dat$sample[grep("RD-20180817-001-SR18271",sub.dat$sample)] <- "lesion1"
sub.dat$sample[grep("RD-20180817-002-SR18271",sub.dat$sample)] <- "lesion2"
sub.dat$sample[grep("T-Bsc1",sub.dat$sample)] <- "T_Bsc1"
#===============================================================================
# run cellphoneDB for each sample 
#===============================================================================
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

sample_name <- unique(sub.dat$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/",sample_name[i],"/")
    tmp <- subset(sub.dat,cells=which(sub.dat$sample==sample_name[i]))
    cellphoneDB_input(tmp,"type.refine",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

#==================================================================
# run in shell 
for i in `ls /public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/MMyeloid/CCC_Data/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done
#===================================================================


#===================================================================
# get analysis result 
# malignant to Myeloid 
# cellphoneDB_Data use Malignant BMDM MG and T cell to run 
#===================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
Malignant.BMDM <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","malignant.BMDM")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".malignant.BMDM") 
	Malignant.BMDM[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Malignant.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Malignant.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) # get co_id for intersection pair id 

Malignant.BMDM.res <- cbind(Malignant.BMDM[[1]][which(Malignant.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Malignant.BMDM[[2]][which(Malignant.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Malignant.BMDM)){
	Malignant.BMDM.res <- cbind(Malignant.BMDM.res,Malignant.BMDM[[i]][which(Malignant.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.BMDM.res <- Malignant.BMDM.res[,c(1,grep("BMDM",colnames(Malignant.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.BMDM.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
Com.res.f <- Com.res.f[,c("A20190305.malignant.BMDM","A20190312.malignant.BMDM","T_Bsc1.malignant.BMDM","lesion1.malignant.BMDM","lesion2.malignant.BMDM",
	"BT1296.malignant.BMDM","BT1297.malignant.BMDM","scrBT1431m.malignant.BMDM","scrBT1432m.malignant.BMDM")]
# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.pvalue.f.RDS")
#=====================================================================================================================
pdf("/public/workspace/lily/Lung2Brain/TME/CCC_Malignant.BMDM.pdf",height=16)
pheatmap::pheatmap(Com.res.f,cluster_cols=F)
dev.off()



#=====================================================================================================================
# means 
#=====================================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
Malignant.BMDM <- list()
for(i in 1:length(sample_name)){
	#tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.mean.f <- tmp.mean[,c("id_cp_interaction","malignant.BMDM")]
	colnames(tmp.mean.f)[2] <- paste0(sample_name[i],".malignant.BMDM") 
	Malignant.BMDM[[i]] <- tmp.mean.f
}
#===========================================================================
tmp.id <- c()
for(i in 1:length(Malignant.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Malignant.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) 
#===============================================
Malignant.BMDM.res <- cbind(Malignant.BMDM[[1]][which(Malignant.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Malignant.BMDM[[2]][which(Malignant.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(Malignant.BMDM)){
	Malignant.BMDM.res <- cbind(Malignant.BMDM.res,Malignant.BMDM[[i]][which(Malignant.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.BMDM.res <- Malignant.BMDM.res[,c(1,grep("BMDM",colnames(Malignant.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.BMDM.res,tmp.mean[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res.f <- Com.res[,c("A20190305.malignant.BMDM","A20190312.malignant.BMDM","T_Bsc1.malignant.BMDM","lesion1.malignant.BMDM","lesion2.malignant.BMDM",
	"BT1296.malignant.BMDM","BT1297.malignant.BMDM","scrBT1431m.malignant.BMDM","scrBT1432m.malignant.BMDM")]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.mean.f.RDS")


#==============================================================================================================================================
# use pvalue and means to melt into 

pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.mean.f.RDS")
# means.res <- means.res[-which(rowSums(means.res)<9),] #this filter is also important ??

# should make a id to merge different samples 
library(reshape)
pvalue <- melt(as.matrix(pvalue.res),id="col.names")
colnames(pvalue) <- c("gene_pair","samples","pvalue")
pvalue$id <- paste0(pvalue$gene_pair,".",pvalue$samples)
# should make a id to merge different samples
mean <- melt(as.matrix(means.res),id="col.names")
colnames(mean) <- c("gene_pair","samples","means_value")
mean$id <- paste0(mean$gene_pair,".",mean$samples)
#========================================================================
merge(mean,pvalue,by="id") -> res.final
res.final$Samples <- sapply(as.vector(res.final$samples.x),function(x){strsplit(x,"\\.")[[1]][1]})
res.final$Samples <- factor(res.final$Samples,levels=c("lesion1","lesion2","A20190305","A20190312","T_Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))

# lly think here to filter is more good 
# just use plot 
#========================================================================
#res.final.f <- res.final[which(res.final$means_value>1),]
res.final.f <- res.final
res.final.f$means_value[which(res.final.f$pvalue==0)] <- 0


pdf("/public/workspace/lily/Lung2Brain/TME/result_plot/Malignant.BMDM_bubble.pdf",useDingbats=F)

library(ggplot2)
ggplot(res.final.f, aes(x=Samples, y=gene_pair.y, size=pvalue)) + geom_point(aes(colour = means_value))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	#scale_colour_gradientn(low = "ForestGreen",mid="white",high = "#F29403")+
	#scale_colour_gradientn(colors=c('white','#F29403'),values=c(0,0.5,20))+
	scale_colour_gradientn(colours = colorRampPalette(c('white','#F29403'))(500))+
	theme(panel.grid.major = element_blank())

dev.off()

#========================================================
# check result 
# want to know why select CXCL2 and ANXA1 
# 2020-12-1 lly think this is more good
# top 20
#========================================================

pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.BMDM.mean.f.RDS")
# set some filter 
# means.res <- means.res[-which(rowSums(means.res)<9),] #this filter is also important
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]
# tmp <- apply(means.res,1,function(x){length(which(x[]>0.5))})



# should make a id to merge different samples 
library(reshape)
pvalue <- melt(as.matrix(pvalue.res),id="col.names")
colnames(pvalue) <- c("gene_pair","samples","pvalue")
pvalue$id <- paste0(pvalue$gene_pair,".",pvalue$samples)
# should make a id to merge different samples
mean <- melt(as.matrix(means.res),id="col.names")
colnames(mean) <- c("gene_pair","samples","means_value")
mean$id <- paste0(mean$gene_pair,".",mean$samples)
#========================================================================
merge(mean,pvalue,by="id") -> res.final
res.final$Samples <- sapply(as.vector(res.final$samples.x),function(x){strsplit(x,"\\.")[[1]][1]})
res.final$Samples <- factor(res.final$Samples,levels=c("lesion1","lesion2","A20190305","A20190312","T_Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))
res.final.f <- res.final
res.final.f$means_value[which(res.final.f$pvalue==0)] <- 0 #pvalue < 0.05  means to 0
res.final.f$type <- "unsure"
res.final.f$type[which(res.final.f$Samples%in%c("A20190305","A20190312","T_Bsc1"))] <- "LCBM"
res.final.f$type[which(res.final.f$Samples%in%c("BT1296","BT1297","scrBT1431m","scrBT1432m"))] <- "LC"
res.final.f$type[which(res.final.f$Samples%in%c("lesion1","lesion2"))] <- "GBM"
#========================================================================
pairs <- as.vector(unique(res.final.f$gene_pair.y))
res.mat <- matrix(NA,ncol=5)
for(i in 1:length(pairs)){
	tmp.pair <- pairs[i]
	tmp1 <- mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LC"&res.final.f$gene_pair.y==tmp.pair)])
	tmp2 <- mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="GBM"&res.final.f$gene_pair.y==tmp.pair)])
	tmp3 <- mean(res.final.f$means_value[which(res.final.f$type=="LC"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])
	tmp4 <- mean(res.final.f$means_value[which(res.final.f$type=="GBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])
	res.mat <- rbind(res.mat,c(tmp.pair,as.numeric(tmp1),as.numeric(tmp2),as.numeric(tmp3),as.numeric(tmp4)))
}

res.mat <- as.data.frame(res.mat)
res.mat <- res.mat[-1,]
colnames(res.mat) <- c("pair","LCBM.LC","LCBM.GBM","LC.LCBM","GBM.LCBM")
res.data <- data.frame(apply(res.mat[,2:5],2,function(x){as.numeric(as.vector(x))}))
res.data$pair <- res.mat$pair
res.data <- res.data[order(res.data$LCBM.LC,decreasing=T),]

#head(res.data[order(res.data$LCBM.LC,decreasing=T),])
#=============================================================================================================
# plot result 
#=============================================================================================================
res.final.ff  <- res.final.f[which(res.final.f$gene_pair.y%in%res.data$pair[1:20]),]
res.final.ff$gene_pair.y <- factor(res.final.ff$gene_pair.y,levels=res.data$pair[1:20])

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Tumor2BMDM_CCC.pdf",useDingbats=F)
library(ggplot2)
ggplot(res.final.ff, aes(x=Samples, y=gene_pair.y, size=pvalue)) + geom_point(aes(colour = means_value))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	#scale_colour_gradientn(low = "ForestGreen",mid="white",high = "#F29403")+
	#scale_colour_gradientn(colors=c('white','#F29403'),values=c(0,0.5,20))+
	scale_colour_gradientn(colours = colorRampPalette(c('white','#F29403'))(500))+
	theme(panel.grid.major = element_blank())

library(pheatmap)
rownames(res.data) <- res.data$pair
pheatmap(res.data[1:20,1,drop=F],cluster_cols=F,cluster_rows=F)
dev.off()



#=============================================================================================================
# use TCGA LUAD to verfiy ANXA1 and CXCL2  with BMS scoore 
#=============================================================================================================
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")

luad_mod <- data.frame(luad_mod)
luad_mod$ANXA1 <- as.numeric(dat["ANXA1",])
luad_mod$CXCL2 <- as.numeric(dat["CXCL2",])


library(ggExtra)
library(ggplot2)
library(ggpubr)

    ylab("ANXA1 expression") + xlab('score') + stat_smooth(method="lm",se=T) + 2
 0   stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()

#=============================================================================================================
# CXCL2 
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CXCL2_TCGA_BMS_cor.pdf",useDingbats=F)
p<-ggplot(luad_mod,aes(x=BMS_test_norm,y=CXCL2)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("CXCL2 expression") + xlab('score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 

dev.off()


# verfiy CCL4
#=============================================================================================================
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")

luad_mod <- data.frame(luad_mod)
luad_mod$CCL4 <- as.numeric(dat["CCL4",])


library(ggExtra)
library(ggplot2)
library(ggpubr)

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCL4_TCGA_BMS_cor.pdf",useDingbats=F)
p<-ggplot(luad_mod,aes(x=BMS_test_norm,y=CCL4)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("CCL4 expression") + xlab('score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()






#======================================================================================================
# other data to verify CCL4 
# 2020-12-8
#======================================================================================================
dat <- read.csv("~/metastasis/data/verify/GSE137762/GSE137762_20190910_allSamples_SchulzEtal_MYTMEinBRM.csv",sep=";")
# transorm to expression data 
info <- read.table("/public/workspace/lily/REF/mmGRCm38.genelength.txt",sep="\t",header=T)
colnames(info) <- c("chr","type","start","end","gene_id","gene_name")
res <- merge(info[,3:6],dat,by.x="gene_id",by.y="ensemblgene_ID")
res$length <- res$end -res$start
# to TPM 
tmp <- t(apply(res[,6:114],1,function(x){x/x[length(x)]}))*10^3
tmp.res <- apply(tmp,2,function(x){x/sum(x)})*10^6
tmp.res <- data.frame(tmp.res)
tmp.res$gene_name <- res$mouse_gene_name

aggregate(.~gene_name,data=tmp.res,FUN=median) -> res.final
rownames(res.final) <- res.final$gene_name
res.final$gene_name <- NULL
res.final$length <- NULL
res.final <- as.matrix(res.final)
saveRDS(res.final,file="/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
# check CCL4 
#====================================================================================================
ccl4.res <- as.data.frame(res.final["Ccl4",])
ccl4.res$type <- "unknow"
colnames(ccl4.res)[1] <- "exp"
ccl4.res$type <- gsub("[0-9]$","",rownames(ccl4.res))

#=====================================================================================================
# plot result 
#=====================================================================================================
# plot ccl4 result 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
ccl4.res <- as.data.frame(dat["Ccl4",])
ccl4.res$type <- "unknow"
colnames(ccl4.res)[1] <- "exp"
ccl4.res$type <- gsub("[0-9]$","",rownames(ccl4.res))

#====================================================================================================
tmp.f <- ccl4.res[grep("BMDM|Blood",rownames(ccl4.res)),]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM","X5x2_D3_BMDM","X5x2_D5_BMDM","X5x2_D10_BMDM","X1x10_D3_BMDM"))


compare_list <- list(c("Control_BloodMonocyte","small_BMDM"),
	c("Control_BloodMonocyte","large_BMDM"))
# ggplot
library(ggplot2)
library(ggpubr)

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCL4_Mouse_varify.pdf",useDingbats=F,width=10)
ggplot(tmp.f,aes(x=type,y=exp,fill=type)) + geom_boxplot(outlier.shape=NA) +  stat_compare_means(comparisons = compare_list)
dev.off()




#==================================================================================================
# scRNA 
# show no difference 
#==================================================================================================
options(stringsAsFactors=F)
dat776 <- read.delim2("/public/workspace/lily/metastasis/data/verify/GSE137762/GSM4080777_20190826_GSH-MS-002TRC.csv",sep=";")
rownames(dat776) <- sapply(strsplit(as.vector(dat776$GENEID),"__"),function(x){x[1]})
dat776$GENEID <- NULL
colnames(dat776) <- paste0("Cell777_",1:ncol(dat776))

# make a seurat object 
library(Seurat)
Obj776 <- CreateSeuratObject(
counts=as.matrix(dat776),
project = "Sample777",
assay = "RNA",
min.cells = 0,
min.features = 0,
names.field = 1,
names.delim = "_",
meta.data = NULL
)
tmp_dat <- Obj776
tmp_dat = NormalizeData(object = tmp_dat)
tmp_dat <- FindVariableFeatures(object = tmp_dat)
# scaling
all.genes <- rownames(x = tmp_dat)
tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
# PCA
tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
# clustering
tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
# T-SNE
tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
saveRDS(tmp_dat,file="/public/workspace/lily/metastasis/data/verify/GSE137762/GSM776.RDS")
saveRDS(tmp_dat,file="/public/workspace/lily/metastasis/data/verify/GSE137762/GSM777.RDS") # the same prepare as 776 sample 

#==========================================================================================
# integrated 
dat1 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSM776.RDS")
dat2 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSM777.RDS")
integration.anchors <- FindIntegrationAnchors(object.list = c(dat1,dat2))
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



#=================================================================================
# try to use nichenet find CCL4 pathway traget gene 
# 2020-12-11
#=================================================================================
dat <- readRDS("/public/workspace/lily/software/nichenet/ligand_target_matrix.rds")





















#==============================================================================================================
# immunecell AI 
# for LUAD 
#==============================================================================================================
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")
dat <- rbind(dat,rep("Group1",ncol(dat)))
dat <- dat[c("20531",rownames(dat)[-20351]),]
rownames(dat)[1] <-  "group"
write.table(dat,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_imuncellAI.txt",sep="\t",row.names=T,col.names=T,quote=F)






#=============================================================================================================
# run CSOmap 
#=============================================================================================================
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo,lib.loc="/public/workspace/lijie/R/x86_64-pc-linux-gnu-library/4.0.2")

#===========================
# load function
# #=============================================================================================================
# dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
# tmp.dat <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type%in%c("maliganant","Myeloid")))
# Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
# tmp.mye <- subset(Myeloid,cells=which(Myeloid$type_group=="LCBM"))

# tmp.dat$CSOmap.type <- "unknow"
# tmp.dat$CSOmap.type[which(tmp.dat$type=="maliganant")] <- "malignant"
# tmp.dat$CSOmap.type[which(colnames(tmp.dat)%in%gsub("_1$","",colnames(tmp.mye)))] <- tmp.mye$type.refine
# tmp.res <- subset(tmp.dat,cells=which(tmp.dat$CSOmap.type%in%c("malignant","BMDM")))
# saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/LCBM_BMDM_tumor.RDS")

# make all cell type 
# another try 
# maybe this way is a better ways 
#==============================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tmp.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
tmp.mye <- subset(Myeloid,cells=which(Myeloid$type_group=="LCBM"))

tmp.dat$CSOmap.type <- tmp.dat$type
tmp.dat$CSOmap.type[which(tmp.dat$type=="maliganant")] <- "malignant"
tmp.dat$CSOmap.type[which(colnames(tmp.dat)%in%gsub("_1$","",colnames(tmp.mye)))] <- tmp.mye$type.refine
saveRDS(tmp.dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/LCBM_AllCell_CSOmap.RDS")

# sampling 
set.seed(12345)
tmp.res <- subset(tmp.dat,cells=colnames(tmp.dat)[sample(1:ncol(tmp.dat),size=1000)])
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/LCBM_cell_1000.RDS")

# table(tmp.res$CSOmap.type)

#      BMDM malignant
#      7433      8243





#============================================================================================================
runCSOmap <- function (seurat.obj, labels, output.file=NULL) {
  labelData <- data.frame(cells=colnames(seurat.obj),
                          labels=as.vector(seurat.obj[[labels, drop=T]]))
  rownames(labelData) <- labelData$cells

  TPM <- as.matrix(seurat.obj$RNA@data)

  #Calculate optimized 3D coordinates
  affinityMat <- getAffinityMat(TPM, LR, verbose = T)

  coords_res <- runExactTSNE_R(
    X = affinityMat,
    no_dims = 3,
    max_iter = 1000,
    verbose = T
  )
  coords <- coords_res$Y
  rownames(coords) <- colnames(TPM)
  colnames(coords) <- c('x', 'y', 'z')

  #Visualization
  #arrange data
  coords_tbl <- bind_cols(cellName = rownames(coords), as.data.frame(coords))

  join_vec <- setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
  cellinfo_tbl <- left_join(coords_tbl, labelData, by = join_vec)

  density_obj <- getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
  cellinfo_tbl <- cellinfo_tbl %>% mutate(density = density_obj)

  #p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")

  #Get significance
  #signif_results <- getSignificance(coords, labels = cellinfo_tbl$labels, verbose = F)
  #contribution_list <- getContribution(TPM, LR, signif_results$detailed_connections)

  res <- list(coords,cellinfo_tbl)
  saveRDS(res,file=paste0(output.file))
}


#=============================================================================================================
#  run_CSOmapR <- function(seuratObj, label, outfile) {
library(dplyr)
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo,lib.loc="/public/workspace/lijie/R/x86_64-pc-linux-gnu-library/4.0.2")

runCSOmap(subset.dat,"CSOmap.type",output.file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/sample1000.RDS")

runCSOmap(dat,"CSOmap.type",output.file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_malignant_CSOmap.RDS")

# 2020-12-2 run all cell type and  all cells 
runCSOmap(dat,"CSOmap.type",output.file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_ALL_CELl_CSOmap.RDS")

#=============================================================================================================
# vislize result 
# 2020-12-2
#=============================================================================================================
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/LCBM_BMDM_tumor.RDS")
# dat$CSOmap.info  <- dat$CSOmap.type
# data <- dat[['RNA']]@data
# # dat$ANXA1[which(data["ANXA1",]>0 & dat$CSOmap.type=="malignant")] <- "malignant_ANXA1+" # almost all malignant cell express ANXA1 gene 
# dat$CSOmap.info[which(data["FPR1",]>0 & dat$CSOmap.type=="BMDM")] <- "BMDM_FPR1+"
# dat$CSOmap.info[which(data["FPR3",]>0 & dat$CSOmap.type=="BMDM")] <- "BMDM_FPR3+"
# # CXCL2 
# dat$CSOmap.info2 <- dat$CSOmap.type
# dat$CSOmap.info2[which(data["CXCL2",]>0 & dat$CSOmap.type=="malignant")] <- "Malignant_CXCL2+"
# dat$CSOmap.info2[which(data["DPP4",]>0 & dat$CSOmap.type=="BMDM")] <- "BMDM_DPP4+"
# # save info 
# tmp <- as.data.frame(dat@meta.data[,c("CSOmap.info","CSOmap.info2")])
# saveRDS(tmp,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_malignant_CSOmap_info.RDS")

#=======================================================================================================
# run in local 
#=======================================================================================================
dat <- readRDS("F:/lly/Lung2Brain/Fig/version_11_16/Myeloid/BMDM_malignant_CSOmap.RDS")
info <- readRDS("F:/lly/Lung2Brain/Fig/version_11_16/Myeloid/BMDM_malignant_CSOmap_info.RDS")
coords <- dat[[1]]
cellinfo_tbl <- dat[[2]]
cellinfo_tbl$CSOmap.info <- info$CSOmap.info
cellinfo_tbl$CSOmap.info2 <- info$CSOmap.info2


# one way to plot 
library(plotly)
set.seed(12345)
plot_ly(cellinfo_tbl, x = ~x, y = ~y, z = ~z, color = ~CSOmap.info2 , text = ~CSOmap.info,colors = c("#c9c9c9","#2db928","#fcd000","#d13814"),size=0.3)


# BMDM       BMDM_FPR1+       BMDM_FPR3+        malignant
#             2010             2410             3013              241
# malignant_ANXA1+
#             8002


#=======================================================================================================
# maybe this plot is more good 
# 2020-12-2
# another way to plot 
#=======================================================================================================
cex=0.8
phi=130
theta=-80
stat = aggregate(cbind(x,y,z)~labels, cellinfo_tbl,mean)
rownames(stat) = stat[,1]
stat=stat[,-1]
dist(stat,)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='BMDM'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#f2af00",
    xlim=range(cellinfo_tbl $x),ylim=range(cellinfo_tbl $y),
    zlim=range(cellinfo_tbl $z))   
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='malignant'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#ce1126",add=T)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='MG'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#5482ab",add=T)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='T_cell'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#7ab800",add=T)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels%in%c('B_cell',"Endothelial","Fibroblast","Oligodendrocyte","unknow")),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#eeeeee",add=T)






#=================================================================
# calculate distance between BMDM and Malignant 
# 2020-12-7
#=================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_ALL_CELL_CSOmap.RDS")
coords <- dat[[1]]
cellinfo_tbl <- dat[[2]]
tmp <- aggregate(.~labels,data=cellinfo_tbl[,2:5],FUN=median)

rownames(tmp) <- tmp$labels
tmp$labels <- NULL
as.matrix(dist(tmp)) -> res


# res plot 
#=================================================================
res[,"malignant"] -> a
a[-5] -> b
barplot(-log2(b))
-log2(b)
#          B_cell            BMDM     Endothelial      Fibroblast              MG
#       4.0418951       4.9311667       6.2920510       6.1086057       4.2232215
# Oligodendrocyte          T_cell          unknow
#       4.8778778       3.6933448       0.8616857

#=================================================================
#
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/distance_Tumor2type.pdf")
barplot(c(4.9311667,4.436423,4.2232215,3.6933448),ylim=c(2,5))
dev.off()

tmp.res <- c(0.03277713,0.05354065,0.07730230,median(c(0.06071113,0.01276156,0.01449194,0.55030917)))
res.f <- as.vector(scale(tmp.res,center=F))
names(res.f) <- c("BMDM","MG","T cell","Others")
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/distance_Tumor2type_new.pdf")
barplot(res.f,ylim=c(0,1.5))
dev.off()

#========================================================================================================================================================
# check Malignant to MG communication  
# 2020-11-30
#========================================================================================================================================================
# pvalue calculate result 
#========================================================================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
Malignant.MG <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","malignant.MG")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".malignant.MG") 
	Malignant.MG[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Malignant.MG)){
	tmp.id <- c(tmp.id,as.vector(Malignant.MG[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) # get co_id for intersection pair id 

Malignant.MG.res <- cbind(Malignant.MG[[1]][which(Malignant.MG[[1]]$id_cp_interaction%in%co_id),],
							Malignant.MG[[2]][which(Malignant.MG[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Malignant.MG)){
	Malignant.MG.res <- cbind(Malignant.MG.res,Malignant.MG[[i]][which(Malignant.MG[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.MG.res <- Malignant.MG.res[,c(1,grep("MG",colnames(Malignant.MG.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.MG.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("MG",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
Com.res.f <- Com.res.f[,c("A20190305.malignant.MG","A20190312.malignant.MG","T_Bsc1.malignant.MG","lesion1.malignant.MG","lesion2.malignant.MG",
	"BT1296.malignant.MG","BT1297.malignant.MG","scrBT1431m.malignant.MG","scrBT1432m.malignant.MG")]
# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Malignant.MG.pvalue.f.RDS")
#=====================================================================================================================
# pdf("/public/workspace/lily/Lung2Brain/TME/CCC_Malignant.MG.pdf",height=16)
# pheatmap::pheatmap(Com.res.f,cluster_cols=F)
# dev.off()



#=====================================================================================================================
# means 
#=====================================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
Malignant.MG <- list()
for(i in 1:length(sample_name)){
	#tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.mean.f <- tmp.mean[,c("id_cp_interaction","malignant.MG")]
	colnames(tmp.mean.f)[2] <- paste0(sample_name[i],".malignant.MG") 
	Malignant.MG[[i]] <- tmp.mean.f
}
#===========================================================================
tmp.id <- c()
for(i in 1:length(Malignant.MG)){
	tmp.id <- c(tmp.id,as.vector(Malignant.MG[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) 
#===============================================
Malignant.MG.res <- cbind(Malignant.MG[[1]][which(Malignant.MG[[1]]$id_cp_interaction%in%co_id),],
							Malignant.MG[[2]][which(Malignant.MG[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(Malignant.MG)){
	Malignant.MG.res <- cbind(Malignant.MG.res,Malignant.MG[[i]][which(Malignant.MG[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.MG.res <- Malignant.MG.res[,c(1,grep("MG",colnames(Malignant.MG.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.MG.res,tmp.mean[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair
Com.res <- Communcation.res[,grep("MG",colnames(Communcation.res))]
Com.res.f <- Com.res[,c("A20190305.malignant.MG","A20190312.malignant.MG","T_Bsc1.malignant.MG","lesion1.malignant.MG","lesion2.malignant.MG",
	"BT1296.malignant.MG","BT1297.malignant.MG","scrBT1431m.malignant.MG","scrBT1432m.malignant.MG")]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Malignant.MG.mean.f.RDS")


#==============================================================================================================================================
# use pvalue and means to melt into 

pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.MG.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Malignant.MG.mean.f.RDS")
means.res <- means.res[-which(rowSums(means.res)<9),]

# should make a id to merge different samples 
library(reshape)
pvalue <- melt(as.matrix(pvalue.res),id="col.names")
colnames(pvalue) <- c("gene_pair","samples","pvalue")
pvalue$id <- paste0(pvalue$gene_pair,".",pvalue$samples)
# should make a id to merge different samples
mean <- melt(as.matrix(means.res),id="col.names")
colnames(mean) <- c("gene_pair","samples","means_value")
mean$id <- paste0(mean$gene_pair,".",mean$samples)
#========================================================================
merge(mean,pvalue,by="id") -> res.final
res.final$Samples <- sapply(as.vector(res.final$samples.x),function(x){strsplit(x,"\\.")[[1]][1]})
res.final$Samples <- factor(res.final$Samples,levels=c("lesion1","lesion2","A20190305","A20190312","T_Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))

library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/TME/result_plot/Malignant.MG_bubble.pdf",useDingbats=F)
ggplot(res.final, aes(x=Samples, y=gene_pair.y, size=pvalue)) + geom_point(aes(colour = means_value))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	#scale_colour_gradientn(low = "ForestGreen",mid="white",high = "#F29403")+
	scale_colour_gradientn(colours = colorRampPalette(c('white','#F29403'))(500))+
	theme(panel.grid.major = element_blank())

dev.off()






#====================================================================================================
# calculate TCGA M1/M2 
# use cibersort M1 M2 marker gene ?
# just use xiaoxuej 
#====================================================================================================
# make a mod 
# tmp <- read.table("/public/workspace/lily/metastasis/data/verify/LM22.txt",header=T,sep="\t")
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("M1_marker","M2_marker"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)
















#===================================================================================================================================
# check BMS to Tumor cellphoneDB result 
# 2020-12-7
#===================================================================================================================================

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
BMDM.Malignant <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","BMDM.malignant")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".BMDM.malignant") 
	BMDM.Malignant[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(BMDM.Malignant)){
	tmp.id <- c(tmp.id,as.vector(BMDM.Malignant[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) # get co_id for intersection pair id 

BMDM.Malignant.res <- cbind(BMDM.Malignant[[1]][which(BMDM.Malignant[[1]]$id_cp_interaction%in%co_id),],
							BMDM.Malignant[[2]][which(BMDM.Malignant[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(BMDM.Malignant)){
	BMDM.Malignant.res <- cbind(BMDM.Malignant.res,BMDM.Malignant[[i]][which(BMDM.Malignant[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
BMDM.Malignant.res <- BMDM.Malignant.res[,c(1,grep("BMDM",colnames(BMDM.Malignant.res)))]
# add pair gene information 
Communcation.res <- merge(BMDM.Malignant.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
Com.res.f <- Com.res.f[,c("A20190305.BMDM.malignant","A20190312.BMDM.malignant","T_Bsc1.BMDM.malignant","lesion1.BMDM.malignant","lesion2.BMDM.malignant",
	"BT1296.BMDM.malignant","BT1297.BMDM.malignant","scrBT1431m.BMDM.malignant","scrBT1432m.BMDM.malignant")]
# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/BMDM.Malignant.pvalue.f.RDS")
#=====================================================================================================================
# pdf("/public/workspace/lily/Lung2Brain/TME/CCC_BMDM.Malignant.pdf",height=16)
# pheatmap::pheatmap(Com.res.f,cluster_cols=F)
# dev.off()



#=====================================================================================================================
# means 
#=====================================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/")
BMDM.Malignant <- list()
for(i in 1:length(sample_name)){
	#tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.mean.f <- tmp.mean[,c("id_cp_interaction","BMDM.malignant")]
	colnames(tmp.mean.f)[2] <- paste0(sample_name[i],".BMDM.malignant") 
	BMDM.Malignant[[i]] <- tmp.mean.f
}
#===========================================================================
tmp.id <- c()
for(i in 1:length(BMDM.Malignant)){
	tmp.id <- c(tmp.id,as.vector(BMDM.Malignant[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) 
#===============================================
BMDM.Malignant.res <- cbind(BMDM.Malignant[[1]][which(BMDM.Malignant[[1]]$id_cp_interaction%in%co_id),],
							BMDM.Malignant[[2]][which(BMDM.Malignant[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(BMDM.Malignant)){
	BMDM.Malignant.res <- cbind(BMDM.Malignant.res,BMDM.Malignant[[i]][which(BMDM.Malignant[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
BMDM.Malignant.res <- BMDM.Malignant.res[,c(1,grep("BMDM",colnames(BMDM.Malignant.res)))]
# add pair gene information 
Communcation.res <- merge(BMDM.Malignant.res,tmp.mean[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res.f <- Com.res[,c("A20190305.BMDM.malignant","A20190312.BMDM.malignant","T_Bsc1.BMDM.malignant","lesion1.BMDM.malignant","lesion2.BMDM.malignant",
	"BT1296.BMDM.malignant","BT1297.BMDM.malignant","scrBT1431m.BMDM.malignant","scrBT1432m.BMDM.malignant")]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/BMDM.Malignant.mean.f.RDS")



#==================================================================================================================
# result plot 
# 
#==================================================================================================================

pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/BMDM.Malignant.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/BMDM.Malignant.mean.f.RDS")
# set some filter 
# means.res <- means.res[-which(rowSums(means.res)<9),] #this filter is also important
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]
# tmp <- apply(means.res,1,function(x){length(which(x[]>0.5))})



# should make a id to merge different samples 
library(reshape)
pvalue <- melt(as.matrix(pvalue.res),id="col.names")
colnames(pvalue) <- c("gene_pair","samples","pvalue")
pvalue$id <- paste0(pvalue$gene_pair,".",pvalue$samples)
# should make a id to merge different samples
mean <- melt(as.matrix(means.res),id="col.names")
colnames(mean) <- c("gene_pair","samples","means_value")
mean$id <- paste0(mean$gene_pair,".",mean$samples)
#========================================================================
merge(mean,pvalue,by="id") -> res.final
res.final$Samples <- sapply(as.vector(res.final$samples.x),function(x){strsplit(x,"\\.")[[1]][1]})
res.final$Samples <- factor(res.final$Samples,levels=c("lesion1","lesion2","A20190305","A20190312","T_Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))
res.final.f <- res.final
res.final.f$means_value[which(res.final.f$pvalue==0)] <- 0 #pvalue < 0.05  means to 0
res.final.f$type <- "unsure"
res.final.f$type[which(res.final.f$Samples%in%c("A20190305","A20190312","T_Bsc1"))] <- "LCBM"
res.final.f$type[which(res.final.f$Samples%in%c("BT1296","BT1297","scrBT1431m","scrBT1432m"))] <- "LC"
res.final.f$type[which(res.final.f$Samples%in%c("lesion1","lesion2"))] <- "GBM"
#========================================================================
# add another filter to show CCL4 
# filter gene pairs which not show in 3 LCBM samples 
#========================================================================
tmp.BM <- res.final.f[which(res.final.f$type=='LCBM'&res.final.f$means_value>0),]
tmp.LC <- res.final.f[which(res.final.f$type=='LC'&res.final.f$means_value>0),]
res.final.f <- res.final.f[which(res.final.f$gene_pair.y%in%names(which(table(tmp.BM$gene_pair.y)==3))),]
res.final.f <- res.final.f[which(res.final.f$gene_pair.y%in%names(which(table(tmp.LC$gene_pair.y)<=2))),]





#========================================================================
pairs <- as.vector(unique(res.final.f$gene_pair.y))
res.mat <- matrix(NA,ncol=5)
for(i in 1:length(pairs)){
	tmp.pair <- pairs[i]
	tmp1 <- mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LC"&res.final.f$gene_pair.y==tmp.pair)])
	tmp2 <- mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="GBM"&res.final.f$gene_pair.y==tmp.pair)])
	tmp3 <- mean(res.final.f$means_value[which(res.final.f$type=="LC"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])
	tmp4 <- mean(res.final.f$means_value[which(res.final.f$type=="GBM"&res.final.f$gene_pair.y==tmp.pair)])-mean(res.final.f$means_value[which(res.final.f$type=="LCBM"&res.final.f$gene_pair.y==tmp.pair)])
	res.mat <- rbind(res.mat,c(tmp.pair,as.numeric(tmp1),as.numeric(tmp2),as.numeric(tmp3),as.numeric(tmp4)))
}

res.mat <- as.data.frame(res.mat)
res.mat <- res.mat[-1,]
colnames(res.mat) <- c("pair","LCBM.LC","LCBM.GBM","LC.LCBM","GBM.LCBM")
res.data <- data.frame(apply(res.mat[,2:5],2,function(x){as.numeric(as.vector(x))}))
res.data$pair <- res.mat$pair
res.data <- res.data[order(res.data$LCBM.LC,decreasing=T),]

#head(res.data[order(res.data$LCBM.LC,decreasing=T),])
#=============================================================================================================
# plot result 
#=============================================================================================================
res.final.ff  <- res.final.f[which(res.final.f$gene_pair.y%in%res.data$pair[1:20]),]
res.final.ff$gene_pair.y <- factor(res.final.ff$gene_pair.y,levels=res.data$pair[1:20])

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM2Tumor_CCC.pdf",useDingbats=F)
library(ggplot2)
ggplot(res.final.ff, aes(x=Samples, y=gene_pair.y, size=pvalue)) + geom_point(aes(colour = means_value))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	#scale_colour_gradientn(low = "ForestGreen",mid="white",high = "#F29403")+
	#scale_colour_gradientn(colors=c('white','#F29403'),values=c(0,0.5,20))+
	scale_colour_gradientn(colours = colorRampPalette(c('white','#F29403'))(500))+
	theme(panel.grid.major = element_blank())

library(pheatmap)
rownames(res.data) <- res.data$pair
pheatmap(res.data[1:20,1,drop=F],cluster_cols=F,cluster_rows=F)
dev.off()





#============================================================================================================
# 2020-12-30
# check result in Co-Tumor with others immune cells communication
#============================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Both_T.BMDM <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Both_T.BMDM")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Both_T.BMDM") 
	Both_T.BMDM[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Both_T.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Both_T.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 

Both_T.BMDM.res <- cbind(Both_T.BMDM[[1]][which(Both_T.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Both_T.BMDM[[2]][which(Both_T.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(Both_T.BMDM)){
	Both_T.BMDM.res <- cbind(Both_T.BMDM.res,Both_T.BMDM[[i]][which(Both_T.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Both_T.BMDM.res <- Both_T.BMDM.res[,c(1,grep("BMDM",colnames(Both_T.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Both_T.BMDM.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
# Communcation.res <- Communcation.res[-394,]  # change duplicate rownames 
rownames(Communcation.res) <- Communcation.res$interacting_pair
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/TF_CCC_res/Both_T.BMDM.pvalue.f.RDS")

# means matrix is same way to prepare
#############################################################################################################
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Both_T.BMDM <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Both_T.BMDM")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Both_T.BMDM") 
	Both_T.BMDM[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Both_T.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Both_T.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 

Both_T.BMDM.res <- cbind(Both_T.BMDM[[1]][which(Both_T.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Both_T.BMDM[[2]][which(Both_T.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(Both_T.BMDM)){
	Both_T.BMDM.res <- cbind(Both_T.BMDM.res,Both_T.BMDM[[i]][which(Both_T.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Both_T.BMDM.res <- Both_T.BMDM.res[,c(1,grep("BMDM",colnames(Both_T.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Both_T.BMDM.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
# Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR")[1],]  # change duplicate rownames 
rownames(Communcation.res) <- Communcation.res$interacting_pair
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/TF_CCC_res/Both_T.BMDM.means.f.RDS")

##################################################################################################################
# result show 
pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/TF_CCC_res/Both_T.BMDM.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/TF_CCC_res/Both_T.BMDM.means.f.RDS")
# set some filter 
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]

fpair <- names(which(apply(pvalue.res,1,function(x){which(length(which(x>0))>=2)})==1)) #pvalue in 2 samples < 0.05
tmp <- data.frame(apply(means.res,1,function(x){mean(x)})) 
tmp.f <- tmp[fpair,,drop=F]
tmp.f <- tmp.f[order(tmp.f[,1]),,drop=F]

















################################################################################################################################
# 2021-3-12
# SLC7A1 gene 

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 

#tumor <- subset(dat,cells=which(dat$maliganant=="tumor"))
tumor <- subset(dat,cells=which(dat$maliganant=="tumor"&dat$type_group=="LCBM")) 
# whether or not to use LCBM data 
tmp.data <- tumor[["RNA"]]@data
tumor$type.TF <- "other"
tumor$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]==0)] <- "CEBPB"
tumor$type.TF[which(tmp.data["CEBPB",]==0&tmp.data["MYBL2",]>0)] <- "MYBL2"
tumor$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]>0)] <- "CE cells"





# 2. plot BMDM CCL4 expression and four type Tumors to plot heatmap
tumor@active.ident <- factor(tumor$type.TF)
AverageExpression(tumor,features="SLC7A1",assays="RNA")

#         CE cells     CEBPB     MYBL2    other
# SLC7A1 0.1711668 0.1129661 0.1182217 0.114194

# BMDM CCL4 : 11.85748


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CEcell_SLC7A1.pdf",useDingbats=F)
VlnPlot(tumor,features="SLC7A1",group.by="type.TF",pt.size=0)
pheatmap::pheatmap(data.frame(row.names=c("CE","CEBPB","MYBL2","others"),SLC7A1=c(0.1711668,0.1129661,0.1182217,0.114194),BMDM.CCL4=rep(11.85748,4)),
	scale="column",cluster_cols=F,cluster_rows=F)
dev.off()



# load Myeloid data 
myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.mye <- subset(myeloid,cells=which(myeloid$type_group%in%c("LC","LCBM")))
sub.mye@active.ident <- factor(sub.mye$type.refine)
AverageExpression(sub.mye,features="CCL4",assays="RNA")

























#===========================================================================================================================================
# 2021-3-13
# use CSOmap result to plot distance in Both Tumor and T cells 
############################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_ALL_CELL_CSOmap.RDS")
cellinfo_tbl <- dat[[2]]

#################################
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM"))
tmp.data <- sub.dat[['RNA']]@data
sub.dat$type.TF <- sub.dat$type
sub.dat$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]==0&sub.dat$type=="maliganant")] <- "CEBPB"
sub.dat$type.TF[which(tmp.data["CEBPB",]==0&tmp.data["MYBL2",]>0&sub.dat$type=="maliganant")] <- "MYBL2"
sub.dat$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]>0&sub.dat$type=="maliganant")] <- "CE cells"

all(colnames(sub.dat)==cellinfo_tbl$cellName)

cellinfo_tbl$label.new <- as.vector(cellinfo_tbl$labels)
cellinfo_tbl$label.new[which(sub.dat$type.TF=="CEBPB")] <- "CEBPB"
cellinfo_tbl$label.new[which(sub.dat$type.TF=="MYBL2")] <- "MYBL2"
cellinfo_tbl$label.new[which(sub.dat$type.TF=="CE cells")] <- "CE cells"

# calculate result 
tmp <- aggregate(.~label.new,data=cellinfo_tbl[,c(2,3,4,7)],FUN=median)
rownames(tmp) <- tmp$label.new
tmp$label.new <- NULL
as.matrix(dist(tmp)) -> res


















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

#saveRDS(inte, file = "/public/workspace/lily/Lung2Brain/Data/inte_rough.rds")
#PCA
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

#=============================================================================================================================================

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

#========================================================================
dat$type <- "unknow"
dat$type[which(dat$seurat_clusters%in%c(1,2,3,18,39,36,9,30,19,34,40,43,15,28,35,38))] <- "T_cell"
dat$type[which(dat$seurat_clusters%in%c(4,6,12,13,16,7,26,25,14))] <- "Myeloid"
dat$type[which(dat$seurat_clusters%in%c(11,23))] <- "B_cell"
dat$type[which(dat$seurat_clusters%in%c(29,27))] <- "Oligodendrocyte"
dat$type[which(dat$seurat_clusters%in%c(22,33))] <- "Fibroblast"
dat$type[which(dat$seurat_clusters%in%c(32))] <- "Endothelial"
dat$type[which(dat$seurat_clusters%in%c(0,5,10,20,21,37,24,17,8,31))] <- "maliganant"


#=========================================================================
dat$maliganant <- "non-tumor"
dat$maliganant[which(dat$type=="maliganant")] <- "tumor"


data <- dat[['RNA']]@data
mean(apply(data[,which(dat$maliganant=="non-tumor")],2,function(x){length(which(x[]>0))}))

mean(apply(data[,which(dat$maliganant=="tumor")],2,function(x){length(which(x[]>0))}))

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")



#===========================================================================================================================================================================
# choose marker 
#===========================================================================================================================================================================
BM <- subset(dat,cells=which(dat$maliganant=="tumor"&dat$type_group=="LCBM"))
data <- BM[['RNA']]@data
sigtumor <- subset(BM,cells=which(data["PTPRC",]==0))
DefaultAssay(sigtumor) <- "integrated"
#FindVariableFeatures
sigtumor <- FindVariableFeatures(sigtumor, selection.method = "vst")

## Scaling the integratedata
all.genes <- rownames(sigtumor)
sigtumor <- ScaleData(sigtumor, features = all.genes)
# PCA
sigtumor <- RunPCA(sigtumor)
# cluster
sigtumor <- FindNeighbors(sigtumor, dims = 1:20)
sigtumor <- FindClusters(sigtumor, resolution = 1)
# TSNE
# If Umap can not use
sigtumor <- RunTSNE(sigtumor,dims = 1:20)
saveRDS(sigtumor,file="/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")



LC <- subset(dat,cells=which(dat$maliganant=="tumor"&dat$type_group=="LC"))
data <- LC[['RNA']]@data
sigtumor <- subset(LC,cells=which(data["PTPRC",]==0))
DefaultAssay(sigtumor) <- "integrated"
#FindVariableFeatures
sigtumor <- FindVariableFeatures(sigtumor, selection.method = "vst")

##Scaling the integratedata
all.genes <- rownames(sigtumor)
sigtumor <- ScaleData(sigtumor, features = all.genes)
#PCA
sigtumor <- RunPCA(sigtumor)
#cluster
sigtumor <- FindNeighbors(sigtumor, dims = 1:20)
sigtumor <- FindClusters(sigtumor, resolution = 0.5)
#TSNE
# if Umap can not use
sigtumor <- RunTSNE(sigtumor,dims = 1:20)
saveRDS(sigtumor,file="/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")


#========================================================================================================================================================
# one way to calculate 
# 2021-1-2
######################################################################################
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")

data1 <- BM[['RNA']]@data
data2 <- LC[['RNA']]@data
genes <- intersect(rownames(data1),rownames(data2))
data <-as.matrix(cbind(data1[genes,],data2[genes,]))

res <- t(apply(data,1,function(x){
	tmp.BM <- x[1:7496]
	tmp.LC <- x[7497:9403]
	c(	length(which(x[]>0))/length(x),
		length(which(tmp.BM>0))/length(tmp.BM), #pct.BM
		length(which(tmp.LC>0))/length(tmp.LC) #pct.LC
	
		)
	}))
res <- data.frame(res)
colnames(res) <- c("pct.all","pct.BM","pct.LC")
res$pct <- res$pct.BM -res$pct.LC

gene.up <- rownames(res)[which(res$pct>0.05)]
gene.down <- rownames(res)[which(res$pct<(-0.05))]


tmp1 <- data1[c(gene.up,gene.down),]
res_BM <- t(apply(tmp1,2,function(x){
	c(
		length(which(names(x)%in%gene.up&x[]>0))/length(gene.up),
		length(which(names(x)%in%gene.down&x[]>0))/length(gene.down)
		)
	}))

res_BM <- data.frame(res_BM)
colnames(res_BM) <- c("pct.up","pct.down")
res_BM$Cluster <- BM$seurat_clusters
res_BM$pct1 <- res_BM$pct.up-res_BM$pct.down
res_BM$pct2 <- res_BM$pct.up/res_BM$pct.down

aggregate(.~Cluster,data=res_BM,FUN=median) -> tmp1
tmp1[order(tmp1$pct1),]
# so choose cluster 6 and 9 
#===================================================================================================================

tmp2 <- data2[c(gene.up,gene.down),]
res_LC <- t(apply(tmp2,2,function(x){
	c(
		length(which(names(x)%in%gene.up&x[]>0))/length(gene.up),
		length(which(names(x)%in%gene.down&x[]>0))/length(gene.down)
		)
	}))

res_LC <- data.frame(res_LC)
colnames(res_LC) <- c("pct.up","pct.down")
res_LC$Cluster <- LC$seurat_clusters
res_LC$pct1 <- res_LC$pct.up-res_LC$pct.down
res_LC$pct2 <- res_LC$pct.up/res_LC$pct.down

aggregate(.~Cluster,data=res_LC,FUN=median) -> tmp2
tmp2[order(tmp2$pct1),]
# so choose cluster 0 and 11 
########################################################




#==================================================================================================================================
# 2020-11-18
# lly make sure use 6,9 
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(6,9)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(0,11)))
tmp <- merge(subset1,subset2) 
# tmp$type_group <- "LC"
# tmp$type_group[which(tmp$orig.ident%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
tmp@active.ident <- as.factor(tmp$type_group)
geneset_tmp <- FindMarkers(tmp,ident.1="LCBM",ident.2="LC",assay="RNA",logfc.threshold = 0.1)
geneset_tmp$pct <- geneset_tmp$pct.1-geneset_tmp$pct.2
geneset_tmp$pct1 <- geneset_tmp$pct.1/geneset_tmp$pct.2


# FindMarkers(BM,ident.1=c(6,9,1,5),ident.2=c(3,7,8,10),assay="RNA",logfc.threshold = 0) -> gene1
# FindMarkers(LC,ident.1=c(6,7,1,2),ident.2=c(0,11,8,10),assay="RNA",logfc.threshold = 0) -> gene2

# gene1$pct <- gene1$pct.1-gene1$pct.2
# gene2$pct <- gene2$pct.1-gene2$pct.2
# gene1$pct1 <- gene1$pct.1/gene1$pct.2
# gene2$pct1 <- gene2$pct.1/gene2$pct.2


#rownames(geneset_tmp[which(geneset_tmp$avg_logFC>1&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct>0.2&geneset_tmp$pct.2<0.2),])

gene.f <- rownames(geneset_tmp[which(geneset_tmp$avg_logFC>1&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct>0.2&geneset_tmp$pct1>2),])
#================================== do some filter 
# 1. filter MT genes
if(length(grep("^MT-",gene.f))>0){
	gene.ff <- gene.f[-grep("^MT-",gene.f)]
}else{
	gene.ff <- gene.f
}

# 2. filter RP genes 
if(length(grep("^RP",gene.ff))>0){
	gene.fff <- gene.ff[-grep("^RP",gene.ff)]
}else{
	gene.fff <- gene.ff
}


# 3. filter LINC 
if(length(grep("^LINC",gene.ff))>0){
	gene.ffff <- gene.fff[-grep("^LINC",gene.ff)]
}else{
	gene.ffff <- gene.fff
}


#====================================================
# another filter filte LM22 genes 
LMgene <- read.table("/public/workspace/lily/Lung2Brain/Opt_sig/LM22.txt",sep="\t",header=T)
if(length(which(gene.fff%in%LMgene$Gene.symbol))>0){
	gene.ffff <- gene.fff[-which(gene.fff%in%LMgene$Gene.symbol)]
}else{
	gene.ffff <- gene.fff
}

#=======================
# save as BMS_test







#geneset_tmp[which(rownames(geneset_tmp)=="MKI67"),]




#==================================================================================================================================
# BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
# LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
# subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(1,6)))
# subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(0,8)))

# data1 <- subset1[['RNA']]@data
# data2 <- subset2[['RNA']]@data
# data <- as.matrix(cbind(data1,data2))

# res <- t(apply(data,1,function(x){
# 	tmp.BM <- x[1:1461]
# 	tmp.LC <- x[1462:ncol(data)]
# 	c(	length(which(x[]>0))/length(x), #pct.all
# 		length(which(tmp.BM>0))/length(tmp.BM), #pct.BM
# 		length(which(tmp.LC>0))/length(tmp.LC), #pct.LC
# 		wilcox.test(tmp.BM,tmp.LC)$p.value,
# 		mean(tmp.BM)/mean(tmp.LC)
# 		)
# 	}))
# res <- data.frame(res)
# colnames(res) <- c("pct.all","pct.BM","pct.LC","p_val","FC")
# res$p_val_adj <- p.adjust(res$p_val)
# res$pct <- res$pct.BM-res$pct.LC
# res$pct1 <- res$pct.BM/res$pct.LC




# res.f <- res[-which(res$pct.all<0.1),]
# # res.f <- res.f[-which(res.f$pct.BM)]
# length(which(res.f$p_val_adj<0.01&res.f$FC>2&res.f$pct>0.5))


# rownames(gene1[which(gene1$avg_logFC>0&gene1$p_val_adj<0.01&gene1$pct1>1.5),]),


# gene.f <- Reduce(intersect,list(
# 					rownames(gene2[which(gene2$avg_logFC>0&gene2$p_val_adj<0.01&gene2$pct1>1.5),]),
#  					rownames(geneset_tmp[which(geneset_tmp$avg_logFC>0&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct1>1.5),])
#  					))




# gene.f <- intersect(rownames(gene1[which(gene1$avg_logFC>0&gene1$p_val_adj<0.01&gene1$pct>0&gene1$pct.2>0),]),
#  					rownames(gene2[which(gene2$avg_logFC>0&gene2$p_val_adj<0.01&gene2$pct>0&gene2$pct.2>0),]))








#===================================================================================================================
#==================================================================================================================================
# BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
# LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
# dat <- merge(BM,LC)
# dat@active.ident <- as.factor(dat$type_group)
# geneset1 <- FindMarkers(dat,ident.1="LCBM",ident.2="LC",assay="RNA")


# source('~/software/ssGSEA/ssgseaMOD.r')
# geneset1$pct <- geneset1$pct.1/geneset1$pct.2
# gene <- rownames(geneset1[which(geneset1$avg_logFC>0&geneset1$p_val_adj<0.01),])
# # logFC is a good parameter to filter genes
# mod.generate(gene,'BM_gene',out='/public/workspace/lily/Lung2Brain/inte7/BM_gene.mod')



# #=============================================================
# # caculate gene for LC



# #=============================================================
# library(Seurat)
# BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
# #load("/public/workspace/lily/Lung2Brain/inte6/LCBM_BM_gene_mod.RData")
# mod <- mod.analyze2(as.matrix(BM[['RNA']]@data),c('BM_gene'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
# save(mod,file="/public/workspace/lily/Lung2Brain/inte7/LCBM_BM_gene_mod.RData")
# mod <- data.frame(mod)
# mod$Cluster <- BM$seurat_clusters
# aggregate(BM_gene_norm~Cluster,data=mod,FUN=median) -> cluster_score1
# cluster_score1$cell_num <- table(BM$seurat_clusters)
# cluster_score1 <- cluster_score1[order(cluster_score1$BM_gene_norm),]





# LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
# #load("/public/workspace/lily/Lung2Brain/inte6/LC_BM_gene_mod.RData")
# mod <- mod.analyze2(as.matrix(LC[['RNA']]@data),c('BM_gene'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
# save(mod,file="/public/workspace/lily/Lung2Brain/inte7/LC_BM_gene_mod.RData")
# mod <- data.frame(mod)
# mod$Cluster <- LC$seurat_clusters
# aggregate(BM_gene_norm~Cluster,data=mod,FUN=median) -> cluster_score2
# cluster_score2$cell_num <- table(LC$seurat_clusters)
# cluster_score2 <- cluster_score2[order(cluster_score2$BM_gene_norm),]


































# 2021-12-27
# this program is used to find signature 
#============================================================================================================================================
# purify Tumor cells by PTPRC 
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
tmp.data <- tmp.dat[["RNA"]]@data
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Tumor"&tmp.data["PTPRC",]==0))


# split into LCBM and nMLUAD 
# and then re-inte to find good cluster
LCBM <- subset(dat,cells=which(dat$type_group=="LCBM"))
nMLUAD <- subset(dat,cells=which(dat$type_group=="nMLUAD"))


#===========================================================================================================================================
# re-cluster for LCBM
inte.list <- list()
samplelist <- unique(LCBM$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(LCBM,cells=which(LCBM$orig.ident==samplelist[i]))
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
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_tumor_pur.RDS")



#==========================================================================================================================================
# re-cluster for nMLUAD
inte.list <- list()
samplelist <- unique(nMLUAD$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(nMLUAD,cells=which(nMLUAD$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=30)
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
inte <- FindClusters(inte,resolution=0.5)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_tumor_pur.RDS")



#==========================================================================================================================================
# 2021-12-28 
# decide to test MLUAD 
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
tmp.data <- tmp.dat[["RNA"]]@data
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Tumor"&tmp.data["PTPRC",]==0))

# and then re-inte to find good cluster
MLUAD <- subset(dat,cells=which(dat$type_group=="MLUAD"))

recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.1)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)

	return(tmp_dat)
}
DefaultAssay(MLUAD) <- "RNA"
MLUAD <- recluster(MLUAD)

saveRDS(MLUAD,file="/public/workspace/lily/Lung2Brain/Version6/Data/MLUAD_tumor_pur.RDS")









#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
# 1. Find LCBM high signature 
# 2. Find nMLUAD high siganture 
# 3. use LCBM high cluster vs nMLUAD high cluster

##############################################################################################################################################
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_tumor_pur.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_tumor_pur.RDS")
numBM <- ncol(BM)
numLC <- ncol(LC)

data1 <- BM[['RNA']]@data
data2 <- LC[['RNA']]@data
data <-as.matrix(cbind(data1,data2))

res <- t(apply(data,1,function(x){
	tmp.BM <- x[1:numBM]
	tmp.LC <- x[(numBM+1):(numBM+numLC)]
	c(	length(which(x[]>0))/length(x),
		length(which(tmp.BM>0))/length(tmp.BM), #pct.BM
		length(which(tmp.LC>0))/length(tmp.LC) #pct.LC
	
		)
	}))
res <- data.frame(res)
colnames(res) <- c("pct.all","pct.BM","pct.LC")
res$pct1 <- res$pct.BM -res$pct.LC # for high in BM 
res$pct2 <- res$pct.LC -res$pct.BM # for high in Lung
res.f <- res[which(res$pct.all > 0.05),] # gene should be expressed in 5% cells


#=========================================================================================================================================
gene.BM <- rownames(res.f)[which(res.f$pct1>0.0)] # rough BM siganture 
gene.LC <- rownames(res.f)[which(res.f$pct2>0.0)] # rough LC signature 


#=========================================================================================================================================
# now calculate gene expressed percentage in clusters
# for BM
tmp1 <- data1[c(gene.BM,gene.LC),]
res_BM <- t(apply(tmp1,2,function(x){
	c(
		length(which(names(x)%in%gene.BM&x[]>0))/length(gene.BM),
		length(which(names(x)%in%gene.LC&x[]>0))/length(gene.LC)
		)
}))

res_BM <- data.frame(res_BM)
colnames(res_BM) <- c("pct.BM","pct.LC")
res_BM$Cluster <- BM$seurat_clusters
res_BM$pct1 <- res_BM$pct.BM-res_BM$pct.LC
res_BM$pct2 <- res_BM$pct.BM/res_BM$pct.LC

aggregate(.~Cluster,data=res_BM,FUN=median) -> tmp.BM
tmp.BM[order(tmp.BM$pct1),]



#============================================================================================================================================
# 
#============================================================================================================================================
# for MLUAD 
MLUAD <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/MLUAD_tumor_pur.RDS")
tmp.mLUAD <- as.matrix(MLUAD[["RNA"]]@data[c(gene.BM,gene.LC),])
res_mLUAD <- t(apply(tmp.mLUAD,2,function(x){
	c(
		length(which(names(x)%in%gene.BM&x[]>0))/length(gene.BM),
		length(which(names(x)%in%gene.LC&x[]>0))/length(gene.LC)
		)
}))

res_mLUAD <- data.frame(res_mLUAD)
colnames(res_mLUAD) <- c("pct.BM","pct.LC")
res_mLUAD$Cluster <- MLUAD$seurat_clusters
res_mLUAD$pct1 <- res_mLUAD$pct.BM-res_mLUAD$pct.LC
res_mLUAD$pct2 <- res_mLUAD$pct.BM/res_mLUAD$pct.LC

aggregate(.~Cluster,data=res_mLUAD,FUN=median) -> tmp.mLUAD
tmp.mLUAD[order(tmp.mLUAD$pct1),]










#============================================================================================================================================
# 
#============================================================================================================================================
# for nMLUAD 
tmp2 <- data2[c(gene.BM,gene.LC),]
res_LC <- t(apply(tmp2,2,function(x){
	c(
		length(which(names(x)%in%gene.BM&x[]>0))/length(gene.BM),
		length(which(names(x)%in%gene.LC&x[]>0))/length(gene.LC)
		)
	}))

res_LC <- data.frame(res_LC)
colnames(res_LC) <- c("pct.BM","pct.LC")
res_LC$Cluster <- LC$seurat_clusters
res_LC$pct1 <- res_LC$pct.BM-res_LC$pct.LC # find high LUAD signature 
res_LC$pct2 <- res_LC$pct.BM/res_LC$pct.LC

aggregate(.~Cluster,data=res_LC,FUN=median) -> tmp.LC
tmp.LC[order(tmp.LC$pct1),]




##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
# now get cluster subset to find result 
library(Seurat)
# calculate  by myself 
BM <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_tumor_pur.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_tumor_pur.RDS")

subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(6,8,10,14)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(2,8,6)))
tmp <- merge(subset1,subset2) 
tmp.data <- as.matrix(tmp[['RNA']]@data)
idx.bm <- unname(which(tmp$type_group=="LCBM"))
idx.lc <- unname(which(tmp$type_group=="nMLUAD"))

tmp.res <- t(apply(tmp.data,1,function(x){
    bm.exp <- mean(expm1(x[idx.bm]))
    lc.exp <- mean(expm1(x[idx.lc]))
    logfc <- log2(bm.exp/lc.exp)
    pct.bm <- length(which(x[idx.bm]>0))/length(idx.bm)
    pct.lc <- length(which(x[idx.lc]>0))/length(idx.lc)
    p <- wilcox.test(expm1(x[idx.bm]),expm1(x[idx.lc]))$p.value

    c(bm.exp,lc.exp,logfc,pct.bm,pct.lc,p)
}))
tmp.res <- data.frame(tmp.res)
colnames(tmp.res) <- c("BM.exp","LC.exp","log2FC","pct.BM","pct.LC","pvalue")
tmp.res$p.adj <- p.adjust(tmp.res$pvalue)

# do some filter 
res.f <- tmp.res[which(tmp.res$p.adj<0.01&tmp.res$pct.BM>0&tmp.res$pct.LC>0),]
res.f$pct1 <- res.f$pct.BM - res.f$pct.LC
res.f$pct2 <- res.f$pct.BM / res.f$pct.LC


gene <- rownames(res.f[which(res.f$BM.exp>1.5&res.f$log2FC>1&res.f$pct.BM>0.7&res.f$pct.LC<0.4&res.f$pct2>2),])

source("/public/workspace/lily/Lung2Brain/Version6/Signature/Signature_test.r")
sigtest(gene=gene,name="gene_test")


# 2022-1-3
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"BMS_V6",out=paste0("/public/workspace/lily/Lung2Brain/Version6/Signature/BMS_V6.mod")) # make a mod file 
saveRDS(gene,file="/public/workspace/lily/Lung2Brain/Version6/Signature/BMS_V6_gene.RDS")








# use function to iteration
for(pct1 in  c(0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75)){
    for(pct2 in c(1,2,3,4,5)){
        tmp.gene <-  rownames(res.f[which(res.f$BM.exp>0.5&res.f$log2FC>2&res.f$pct1>pct1&res.f$pct2>pct2),]) 
        Lung.gene <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Signature/gene_MLUAD_nMLUAD.RDS")
        gene <- intersect(tmp.gene,Lung.gene)
        source("/public/workspace/lily/Lung2Brain/Version6/Signature/Signature_test.r")
        sigtest(gene=gene,name=paste0("gene_pct1_",pct1,"_pct2_",pct2))
    }
}



#============================================================================================================================================
# 2021-12-28 
# try to use MLUAD vs. nMLUAD to find signature to do intersection
#============================================================================================================================================
library(Seurat)
MLUAD <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/MLUAD_tumor_pur.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_tumor_pur.RDS")

subset1 <- subset(MLUAD,cells=which(MLUAD$seurat_clusters%in%c(1,2)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(2,8,6)))
tmp <- merge(subset1,subset2) 
tmp.data <- as.matrix(tmp[['RNA']]@data)
idx.mluad <- unname(which(tmp$type_group=="MLUAD"))
idx.lc <- unname(which(tmp$type_group=="nMLUAD"))

tmp.res <- t(apply(tmp.data,1,function(x){
    mluad.exp <- mean(expm1(x[idx.mluad]))
    lc.exp <- mean(expm1(x[idx.lc]))
    logfc <- log2(mluad.exp/lc.exp)
    pct.mluad <- length(which(x[idx.mluad]>0))/length(idx.mluad)
    pct.lc <- length(which(x[idx.lc]>0))/length(idx.lc)
    p <- wilcox.test(expm1(x[idx.mluad]),expm1(x[idx.lc]))$p.value

    c(mluad.exp,lc.exp,logfc,pct.mluad,pct.lc,p)
}))
tmp.res <- data.frame(tmp.res)
colnames(tmp.res) <- c("MLUAD.exp","LC.exp","log2FC","pct.MLUAD","pct.LC","pvalue")
tmp.res$p.adj <- p.adjust(tmp.res$pvalue)

# do some filter 
res.f <- tmp.res[which(tmp.res$p.adj<0.01&tmp.res$pct.MLUAD>0&tmp.res$pct.LC>0),]
res.f$pct1 <- res.f$pct.MLUAD - res.f$pct.LC
res.f$pct2 <- res.f$pct.MLUAD / res.f$pct.LC


# use Paired Lung cancer data to filter 
Lung.gene <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Signature/gene_MLUAD_nMLUAD.RDS")

tmp.gene <- merge(Lung.gene[,c("MLUAD.exp","LC.exp","log2FC","pct.MLUAD","pct1","pct2")],
    res.f[,c("BM.exp","pct.LC","log2FC","pct.BM","pct1","pct2")],by="row.names")

rownames(tmp.gene) <- tmp.gene$Row.names
tmp.gene$Row.names <- NULL

colnames(tmp.gene)[3] <- "log2FC.MLUAD"
colnames(tmp.gene)[5] <- "pct1.MLUAD"
colnames(tmp.gene)[6] <- "pct2.MLUAD"
colnames(tmp.gene)[9] <- "log2FC.BM"
colnames(tmp.gene)[11] <- "pct1.BM"
colnames(tmp.gene)[12] <- "pct2.BM"

gene <- rownames(
  tmp.gene[which(tmp.gene$MLUAD.exp>0.5&tmp.gene$pct.LC<0.3&tmp.gene$BM.exp>0.5&tmp.gene$log2FC.MLUAD>0.5&tmp.gene$log2FC.BM>1&tmp.gene$pct.MLUAD>0.5&tmp.gene$pct.BM>0.3),]
  
)

# another way to use pair Lung data to filter 
Lung.gene <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Signature/gene_MLUAD_nMLUAD.RDS")
Lung.gene <- Lung.gene[order(Lung.gene$pct1,decreasing=T),]

gene <- intersect(rownames(res.f[which(res.f$BM.exp>1&res.f$log2FC>1&res.f$pct.BM>0.6&res.f$pct.LC<0.4&res.f$pct1>0.4),]),
 rownames(Lung.gene[which(Lung.gene$MLUAD.exp>2&Lung.gene$log2FC>1&Lung.gene$pct.MLUAD>0.5&Lung.gene$pct.LC<0.4&Lung.gene$pct1>0.2),])
)





source("/public/workspace/lily/Lung2Brain/Version6/Signature/Signature_test.r")
sigtest(gene=gene,name="gene6")


write.table(gene,file="~/tmp.txt",sep="\t",quote=F,row.names=F,col.names=F)


#==================================================================================================================================



gene.Lung <- res.f
saveRDS(gene.Lung,file="/public/workspace/lily/Lung2Brain/Version6/Signature/gene_MLUAD_nMLUAD.RDS")






gene <- rownames(res.f[which(res.f$BM.exp>1.5&res.f$log2FC>1&res.f$pct.BM>0.75&res.f$pct.LC<0.5&res.f$pct2>2),])




# rownames(res.f[which(res.f$BM.exp>1&res.f$log2FC>2&res.f$pct.BM>0.7&res.f$pct.LC<0.5&res.f$pct1>0.5),])


source("/public/workspace/lily/Lung2Brain/Version6/Signature/Signature_test.r")
sigtest(gene=gene,name="gene6")












































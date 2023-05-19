
# 2021-4-29 
# calculate BMS signature 
#####################################################################################################################
# code recording 

# 0. get tumor cell from Inte7_ann.RDS
#####################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tumor <- subset(dat,cells=which(dat$type=="maliganant"&dat$type_group=="LCBM"))
tumor$sample <- tumor$orig.ident
# re clustering 
inte.list <- list()
samplelist <- unique(tumor$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(tumor,cells=which(tumor$sample==samplelist[i]))
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
inte <- FindClusters(inte,resolution=0.5)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/inte7/Data/LCBM_tumor.RDS")
#====================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tumor <- subset(dat,cells=which(dat$type=="maliganant"&dat$type_group=="LC"))
tumor$sample <- tumor$orig.ident
# re clustering 
inte.list <- list()
samplelist <- unique(tumor$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(tumor,cells=which(tumor$sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.score= 100,k.filter= 100)
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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/inte7/Data/LC_tumor.RDS")












# 1. purify tumor cells with PTPRC 
#####################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LCBM_tumor.RDS")
tmp.data <- dat[['RNA']]@data
sub.dat <- subset(dat,cells=which(tmp.data["PTPRC",]==0))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/version_4_27/pur_LCBM_tumor.RDS")
rm(list=ls())
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LC_tumor.RDS")
tmp.data <- dat[['RNA']]@data
sub.dat <- subset(dat,cells=which(tmp.data["PTPRC",]==0))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/version_4_27/pur_LC_tumor.RDS")









# 2. re-integrtion and re-clustering tumor cells 
#####################################################################################################################
# re-clustering
library(Seurat)
lcbm <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/pur_LCBM_tumor.RDS")
inte.list <- list()
samplelist <- unique(lcbm$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(lcbm,cells=which(lcbm$sample==samplelist[i]))
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
inte <- FindClusters(inte,resolution=0.5)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
#=====================================================================================================================
luad <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/pur_LC_tumor.RDS")
inte.list <- list()
samplelist <- unique(luad$sample)
for(i in 1:length(samplelist)){
    tmp <- subset(luad,cells=which(luad$sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.score= 50,k.filter= 50)
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
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")
















# 3. 0 calculate pathway to re-define 
####################################################################################################################
dat <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")

pathway<-c("Apoptosis","Autophagy","Chemokine_signaling_pathway","Circadian_rhythm","Cytosolic_DNA-sensing_pathway","Endocytosis",
  "HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_HYPOXIA","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","KEGG_Focal_adhesion","KEGG_JAK-STAT_signaling_pathway","KEGG_MAPK_signaling_pathway",
  "KEGG_NF-kappa_B_signaling_pathway","KEGG_Notch_signaling_pathway","KEGG_p53_signaling_pathway","KEGG_PI3K-Akt_signaling_pathway",
  "KEGG_TGF-beta_signaling_pathway","KEGG_Wnt_signaling_pathway","Toll-like_receptor_signaling_pathway")

source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),pathway,'/public/workspace/zhangtt/Project/PDAC/Mod/',permN=0)
mod <- as.data.frame(mod)
mod.f <- mod[,22:42]
mod.f$Cluster <- paste0("C",dat$seurat_clusters)
tmp.res <- aggregate(.~Cluster,data=mod.f,FUN=median)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL
res.f <- t(tmp.res)
pdf("/public/workspace/lily/Lung2Brain/version_4_27/diedai/LCBM_metastasis_pathwway.pdf")
pheatmap::pheatmap(res.f,scale="row")
dev.off()
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/version_4_27/diedai/LCBM_metastasis_pathway.RDS")

# hallmark 
#==================================================================================================================
rm(list=ls())
dat <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),modlist,'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- data.frame(mod)
mod.f <- mod[,51:100]
mod.f$Cluster <- paste0("C",dat$seurat_clusters)
tmp.res <- aggregate(.~Cluster,data=mod.f,FUN=median)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL
res.f <- t(tmp.res)
pdf("/public/workspace/lily/Lung2Brain/version_4_27/diedai/LCBM_Hallmark.pdf")
pheatmap::pheatmap(res.f,scale="row")
dev.off()
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/version_4_27/diedai/LCBM_hallmark.RDS")






# 3.1 choose clusters
##################################################################################################################
# now calculate clusters
#
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")

data1 <- BM[['RNA']]@data
data2 <- LC[['RNA']]@data
data <-as.matrix(cbind(data1,data2))

res <- t(apply(data,1,function(x){
	tmp.BM <- x[1:7496]
	tmp.LC <- x[7497:9403]
	c(	length(which(x[]>0))/length(x),
		length(which(tmp.BM>0))/length(tmp.BM), # pct.BM
		length(which(tmp.LC>0))/length(tmp.LC) # pct.LC
	
		)
	}))
res <- data.frame(res)
colnames(res) <- c("pct.all","pct.BM","pct.LC")
res$pct <- res$pct.BM -res$pct.LC
res.f <- res[which(res$pct.all > 0.05),] # pct.all  means this gene need to expression in 5% cells
gene.up <- rownames(res.f)[which(res.f$pct>0)]
gene.down <- rownames(res.f)[which(res.f$pct<0)]

################################################################################################################
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
# so choose cluster 4 and 7 
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
# so choose 1 and 8




# 4. calculate gene 
#######################################################################################################################
# library(Seurat)
# BM <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
# LC <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")
# subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(4,7)))
# subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(1,8)))
# tmp <- merge(subset1,subset2) 
# # tmp$type_group <- "LC"
# # tmp$type_group[which(tmp$orig.ident%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
# tmp@active.ident <- as.factor(tmp$type_group)
# geneset_tmp <- FindMarkers(tmp,ident.1="LCBM",ident.2="LC",assay="RNA",logfc.threshold = 0.1,min.pct = 0.1)
# geneset_tmp$pct <- geneset_tmp$pct.1-geneset_tmp$pct.2
# geneset_tmp$pct1 <- geneset_tmp$pct.1/geneset_tmp$pct.2



# calculate  by myself 
BM <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")
subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(4,7)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(1,8)))
tmp <- merge(subset1,subset2) 
tmp.data <- as.matrix(tmp[['RNA']]@data)
idx.bm <- unname(which(tmp$type_group=="LCBM"))
idx.lc <- unname(which(tmp$type_group=="LC"))

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


res.f <- tmp.res[which(tmp.res$p.adj<0.01&tmp.res$pct.BM>0.01&tmp.res$pct.LC>0.01),]
res.f$pct1 <- res.f$pct.BM - res.f$pct.LC
res.f$pct2 <- res.f$pct.BM / res.f$pct.LC



rownames(res.f[which(res.f$BM.exp>1&res.f$log2FC>2&res.f$pct1>0.25&res.f$pct2>5),]) -> gene


# > gene
#  [1] "F3"       "IER5"     "ASPM"     "G0S2"     "TP53I3"   "BCYRN1"
#  [7] "SFTPB"    "GLS"      "ALCAM"    "MGLL"     "NMU"      "HOPX"
# [13] "IL8"      "CXCL2"    "EREG"     "FAM13A"   "IRX2"     "FCHO2"
# [19] "ELL2"     "TGFBI"    "HMMR"     "MGAT4B"   "SKAP2"    "ANLN"
# [25] "STARD3NL" "UPP1"     "TFPI2"    "BRI3"     "CAV1"     "MET"
# [31] "CREB3L2"  "RARRES2"  "SH3KBP1"  "PLS3"     "CEBPD"    "CHCHD7"
# [37] "LACTB2"   "FAM83A"   "PLEC"     "LCN2"     "PLAU"     "SNCG"
# [43] "MYOF"     "RGS10"    "TAF10"    "SAA1"     "VEGFB"    "BIRC3"
# [49] "SCNN1A"   "GPRC5A"   "LMO3"     "KRT6A"    "DUSP6"    "RGCC"
# [55] "KLF5"     "RNASE1"   "MPRIP"    "ITGA3"    "COL1A1"   "MXRA7"
# [61] "TK1"      "BIRC5"    "BAIAP2"   "ALYREF"   "TPX2"     "TGM2"
# [67] "DMKN"     "C19orf33" "CEACAM6"  "NAPSA"    "COL6A1"   "COL6A2"
# > length(gene)
# [1] 72


gene <- c('F3','IER5','ASPM','G0S2','TP53I3','BCYRN1','SFTPB','GLS','ALCAM','MGLL','NMU','HOPX','IL8','CXCL2','EREG','FAM13A','IRX2','FCHO2','ELL2','TGFBI','HMMR','MGAT4B','SKAP2','ANLN','STARD3NL','UPP1','TFPI2','BRI3','CAV1','MET','CREB3L2','RARRES2','SH3KBP1','PLS3','CEBPD','CHCHD7','LACTB2','FAM83A','PLEC','LCN2','PLAU','SNCG','MYOF','RGS10','TAF10','SAA1','VEGFB','BIRC3','SCNN1A','GPRC5A','LMO3','KRT6A','DUSP6','RGCC','KLF5','RNASE1','MPRIP','ITGA3','COL1A1','MXRA7','TK1','BIRC5','BAIAP2','ALYREF','TPX2','TGM2','DMKN','C19orf33','CEACAM6','NAPSA','COL6A1','COL6A2'
)
saveRDS(gene,file="/public/workspace/lily/Lung2Brain/Version5/Signature/BMS_update_gene.RDS")


# generate a mouse gene
tmp <- read.table("~/REF/Mouse_Human_gene.txt",sep="\t",header=T)
unique(as.vector(tmp[which(tmp$Human.gene.name%in%gene),1])) -> gene.m
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene.m,"BMS_update.m",out=paste0("/public/workspace/lily/Lung2Brain/Version5/Signature/BMS_update.m.mod")) # make a mod file 
mod.generate(gene.m,"BMS_update.m",out=paste0("/public/workspace/lily/MOD_file/BMS_update.m.mod")) # make a mod file

# make a mod file 
#=================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"BMS_update",out=paste0("/public/workspace/lily/Lung2Brain/Version5/Signature/BMS_update.mod")) # make a mod file 

# source("/public/workspace/lily/Lung2Brain/version_4_27/sig_test.r")
# sigtest(gene,name="BMS")


















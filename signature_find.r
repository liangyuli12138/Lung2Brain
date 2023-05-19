
# 2021-4-27
# this program is used to find gene signature 
#########################################################################################################################################
# use hallmark and other result has show C5 maybe the subgroup 

# 0. how LCBM tumor RDS come from . 
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
#=======================================================================================================
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
#===========================================================================================================================================
# try to find gene expression percentage 
# purify tumor cell 
lcbm <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/pur_LCBM_tumor.RDS")
luad <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/pur_LC_tumor.RDS")
index1 <- which(lcbm$seurat_clusters==5) # C5 cells
# index2 <- which(lcbm$seurat_clusters==2) # C2 cells
data1 <- as.matrix(lcbm[['RNA']]@data)
data2 <- as.matrix(luad[["RNA"]]@data)


res <- t(apply(data1,1,function(x){
	tmp.C5 <- x[index1]
	tmp.other <- x[-index1]
	c(	length(which(x[]>0))/length(x),
		length(which(tmp.C5>0))/length(tmp.C5), #pct.C5
		length(which(tmp.other>0))/length(tmp.other) #pct.other
	
		)
	}))
res <- data.frame(res)
colnames(res) <- c("pct.all","pct.C5","pct.other")
res$pct <- res$pct.C5 -res$pct.other
res.bm <- res[order(res$pct,decreasing=T),]
# head(res[which(res$pct.all>0.5&res$pct>0.1),])



# use Lung cancer to filter 
data <-as.matrix(cbind(data1,data2))

res.filter <- t(apply(data,1,function(x){
    tmp.C5 <- x[index1]
	tmp.BM <- x[1:7496]
	tmp.LC <- x[7497:9403]
	c(	length(which(x[]>0))/length(x),
        length(which(tmp.C5>0))/length(tmp.C5),
		length(which(tmp.BM>0))/length(tmp.BM), #pct.BM
		length(which(tmp.LC>0))/length(tmp.LC) #pct.LC

		)
	}))
res.filter <- data.frame(res.filter)
colnames(res.filter) <- c("pct.all","pct.C5","pct.BM","pct.LC")
res.filter$pct1 <- res.filter$pct.C5 -res.filter$pct.LC
res.filter$pct2 <- res.filter$pct.BM -res.filter$pct.LC


# another way to filter 
#============================================================================
# res <- t(apply(data2,1,function(x){
# 	tmp.C2 <- x[index2]
# 	tmp.other <- x[-index2]
# 	c(	length(which(x[]>0))/length(x),
# 		length(which(tmp.C2>0))/length(tmp.C2), #pct.C2
# 		length(which(tmp.other>0))/length(tmp.other) #pct.other
	
# 		)
# 	}))
# res <- data.frame(res)
# colnames(res) <- c("pct.all","pct.C2","pct.other")
# res$pct <- res$pct.C2 -res$pct.other
# res.lc <- res[order(res$pct,decreasing=T),]
# head(res[which(res$pct.all>0.5&res$pct>0.1),])



# for(i in seq(from=0.55,to=0.9,by=0.05)){
#     for(j in seq(from=0.5,to=0.9,by=0.05)){
#         res.f <- res[which(res$pct.all>0.5&res$pct.C5>i&res$pct>0.10),]
#         res.filter.f <- res.filter[which(res.filter$pct1>j&res.filter$pct2>0.5),]
#         gene <- intersect(rownames(res.filter.f),rownames(res.f))

#         if(length(grep("^MT-",gene))>0){
#             gene <- gene[-grep("^MT-",gene)]
#         }
#         if(length(gene)<=1){
#             next
#         }
#         source("/public/workspace/lily/Lung2Brain/version_4_27/sig_test.r")
#         sigtest(gene,name=paste0("Type_0.5_",i,"_",j,"_gene"))

#     }
# }



res.f <- res[which(res$pct.all>0.5&res$pct.C5>0.85&res$pct>0.10),]
res.filter.f <- res.filter[which(res.filter$pct1>0.7&res.filter$pct2>0.5),]
gene <- intersect(rownames(res.filter.f),rownames(res.f))
gene <- gene[-grep("^MT-",gene)]

write.table(gene,file="~/tmp.txt",sep="\t",row.names=F,quote=F,col.names=F)





# test 
source("/public/workspace/lily/Lung2Brain/version_4_27/sig_test.r")
sigtest(gene,name=paste0("Type_0.5_",i,"_",j,"_gene"))









#=====================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LC_tumor.RDS") 
# Calculata Lung cancer hallmark 

pathway<-c("Apoptosis","Autophagy","Chemokine_signaling_pathway","Circadian_rhythm","Cytosolic_DNA-sensing_pathway","Endocytosis",
  "HALLMARK_ANGIOGENESIS","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_HYPOXIA","HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INFLAMMATORY_RESPONSE","KEGG_Focal_adhesion","KEGG_JAK-STAT_signaling_pathway","KEGG_MAPK_signaling_pathway",
  "KEGG_NF-kappa_B_signaling_pathway","KEGG_Notch_signaling_pathway","KEGG_p53_signaling_pathway","KEGG_PI3K-Akt_signaling_pathway",
  "KEGG_TGF-beta_signaling_pathway","KEGG_Wnt_signaling_pathway","Toll-like_receptor_signaling_pathway")

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LC_tumor.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),pathway,'/public/workspace/zhangtt/Project/PDAC/Mod/',permN=0)
mod <- as.data.frame(mod)
mod.f <- mod[,22:42]
mod.f$Cluster <- paste0("C",dat$seurat_clusters)
tmp.res <- aggregate(.~Cluster,data=mod.f,FUN=median)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL
res.f <- t(tmp.res)
pdf("/public/workspace/lily/Lung2Brain/version_4_27/LUAD_metastasis_pathwway.pdf")
pheatmap::pheatmap(res.f,scale="row")
dev.off()

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LC_metastasis_pathway.RDS")





dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LC_tumor.RDS")
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
pdf("/public/workspace/lily/Lung2Brain/version_4_27/LC_Hallmark.pdf")
pheatmap::pheatmap(res.f,scale="row")
dev.off()
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LC_hallmark.RDS")



























#============================================================================================================================================
# diedai re-clustering 
# 2021-4-28
#=============================================================================================================================================
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

##################################################################################################################
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



#=====================================================================================================================
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



# hallmark ========================================================================================================
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
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_hallmark.RDS")























# now calculate clusters
#=====================================================================================================================
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
		length(which(tmp.BM>0))/length(tmp.BM), #pct.BM
		length(which(tmp.LC>0))/length(tmp.LC) #pct.LC
	
		)
	}))
res <- data.frame(res)
colnames(res) <- c("pct.all","pct.BM","pct.LC")
res$pct <- res$pct.BM -res$pct.LC
res.f <- res[which(res$pct.all > 0.05),]
gene.up <- rownames(res.f)[which(res.f$pct>0)]
gene.down <- rownames(res.f)[which(res.f$pct<0)]


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


######################################################################################################################
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")
subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(4,7)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(1,8)))
tmp <- merge(subset1,subset2) 
# tmp$type_group <- "LC"
# tmp$type_group[which(tmp$orig.ident%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
tmp@active.ident <- as.factor(tmp$type_group)
geneset_tmp <- FindMarkers(tmp,ident.1="LCBM",ident.2="LC",assay="RNA",logfc.threshold = 0.1,min.pct = 0.1)
geneset_tmp$pct <- geneset_tmp$pct.1-geneset_tmp$pct.2
geneset_tmp$pct1 <- geneset_tmp$pct.1/geneset_tmp$pct.2



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


gene <- rownames(geneset_tmp[which(geneset_tmp$avg_logFC>1&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct>0.5&geneset_tmp$pct1>1),])

gene <- rownames(geneset_tmp[which(geneset_tmp$avg_logFC>1&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct1>4),])
if(length(grep("^MT-",gene))>0){
	gene <- gene[-grep("^MT-",gene)]
}else{
	gene <- gene
}

write.table(gene,file="~/tmp.txt",sep="\t",row.names=F,quote=F,col.names=F)

# test 
source("/public/workspace/lily/Lung2Brain/version_4_27/sig_test.r")
sigtest(gene,name="BMS")




# test -5
# rownames(res.f[which(res.f$BM.exp>1&res.f$log2FC>2&res.f$pct1>0.25&res.f$pct2>5),]) -> gene







# test for gene_version5
BM <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LCBM_clust.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/version_4_27/diedai/pur_LUAD_clust.RDS")

data.bm <- BM[['RNA']]@data

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(data.bm),c("test_5"),"/public/workspace/lily/Lung2Brain/version_4_27/sig_test/test_5/",permN=0)
mod <- as.data.frame(mod)




data.lc <- LC[["RNA"]]@data
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(data.lc),c("test_5"),"/public/workspace/lily/Lung2Brain/version_4_27/sig_test/test_5/",permN=0)
mod <- as.data.frame(mod)



dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/LCBM_tumor.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[["RNA"]]@data),c("test_5"),"/public/workspace/lily/Lung2Brain/version_4_27/sig_test/test_5/",permN=0)
mod <- as.data.frame(mod)










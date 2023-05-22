

# this program is used to analysis different tumor cell group function
# C0,C9 is Lung gene high 
# C4,C7 is BMS high
# C1,C6 is Brain gene high
# 2021-5-14
# 0. BMS high cells is a group High Plascity cells.
# 1. hallmarks 
# 2. metabolism 
# 3. localinvasion and extravasion 


# 0.0 calculate Brain specific genes # "NEFL" ,"NP22",
gene <- c("GFAP","STMN2","GAP43","PLP1","PMP2",
    "GABBR2","INA","DNM1","SH3GL2","SNAP25",
    "MBP","TPPP","S100B","NEFM","NCAN","KIF3C","SCG3","STMN4")

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$group)

# this show that group16 is like brain 
FeaturePlot(dat,features=c("PLP1","MBP","GAP43","GFAP"),order=T)






###############################################################################################################
# 0. High Plascity cells
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
sub.dat <- subset(dat,cells=which(dat$seurat_clusters%in%c(0,1,4,6,7,9)))

# calculate HPSC signature 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("HPSC_C5"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Data/sub.Tumor.HPSC.RDS")
sub.dat$HPSC <- mod[,2]

saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
# 2021-5-19 
# both Lung gene and BMS is high plascity
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
dat$seurat_clusters <- factor(as.vector(dat$seurat_clusters),levels=c("4","7","0","9","1","6"))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/subtumor_HPSC.pdf",useDingbats=F)
boxplot(HPSC~seurat_clusters,data=dat@meta.data,FUN=median,outline=F)
dev.off()



##################################################################################################################################
# 2021-5-19
# 1. metabolism 
# ~/Lung2Brain/Version5/Metabolism/metabolic_tumor.r.exec
# run sub dat metabolism
# ~/Lung2Brain/Version5/Fig4

# check result 
#=========================================================
# tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Metabolism/KEGGpathway_activity_shuffle.txt",header=T,sep="\t")
# colnames(tmp.dat) <- gsub("X","C",colnames(tmp.dat))

tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/KEGGpathway_activity_shuffle.txt",header=T,sep="\t")
colnames(tmp.dat) <- gsub("X","C",colnames(tmp.dat))
tmp.dat[is.na(tmp.dat)] <- 1

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/subtumor_metabolism.pdf",useDingbats=F,height=10)
pheatmap::pheatmap(tmp.dat,scale="row")
dev.off()





# 2. localvasion and angiogenesis signature 
# looks not every good 
# localinvasion socre calculate 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")

# calculate localinvasion and agogensis
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("localInvasion_cp","intraInvasion_cp","blood_signature"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

tmp.res <- mod[,c(4:6)]
tmp.res$Cluster <- paste0("C",dat$seurat_clusters)




# 3. hallmark calculate 
# 2021-5-20
#============================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat[["RNA"]]@data),modlist,'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/subdat.tumor.hallmark.RDS")

tmp.res <- mod[,51:100]
tmp.res$Cluster <- paste0("C",dat$seurat_clusters)
res.f <- aggregate(.~Cluster,data=tmp.res,FUN=median)
rownames(res.f) <- res.f$Cluster
res.f$Cluster <- NULL
res.f <- t(res.f)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/subtumor_hallmark.pdf",useDingbats=F)
pheatmap::pheatmap(res.f,scale="row")
dev.off()






# 4. MKI67 gene check 
# 2021-5-31
# test if sub tumor 
# 2021-7-1 add other genes
#================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/MKI67_LCBM.pdf",useDingbats=F)
FeaturePlot(dat,feature="MKI67",min.cutoff=0)
FeaturePlot(dat,feature="TOP2A",min.cutoff=0)
FeaturePlot(dat,feature="STMN1",min.cutoff=0)
dev.off()


# use this three gene to calculate a proliferation score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
gene <- c("MKI67","TOP2A","STMN1")
name <- "Proliferation"
respath <- "/public/workspace/lily/MOD_file/NatureMedicine/"
mod.generate(gene,name,out=paste0(respath,name,".mod")) # make a mod file 

# use sub tumor to calculate 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("Proliferation"),"/public/workspace/lily/MOD_file/NatureMedicine/",permN=0)
mod <- as.data.frame(mod)

mod$group <- dat$group
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/porliferation_LCBM_sub.pdf",useDingbats=F)
boxplot(Proliferation_norm~group,data=mod,outline=F,ylim=c(0,1))
dev.off()



# 5. pySCENIC 
# 2021-5-31
#=================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM/step3.auc_mtx.csv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
all(rownames(tf.dat)==rownames(dat@meta.data))

# tf data 
tf.dat$Cluster <- dat$seurat_clusters
tmp.res <- aggregate(.~Cluster,data=tf.dat,FUN=mean)
rownames(tmp.res ) <- paste0("C",tmp.res$Cluster)
tmp.res$Cluster <- NULL

# caculate which TFs 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/TFs_LCBM.pdf",useDingbats=F,width=30,height=10)
pheatmap::pheatmap(tmp.res,scale="column",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()


# try to use C0 C9 C1 C6 C4 C7 
#================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM/step3.auc_mtx.csv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
all(rownames(tf.dat)==rownames(dat@meta.data))

sub.dat <- subset(dat,cells=which(dat$seurat_clusters%in%c(0,9,1,6,4,7)))
sub.dat$group <- "group09"
sub.dat$group[which(sub.dat$seurat_clusters%in%c(1,6))] <- "group16"
sub.dat$group[which(sub.dat$seurat_clusters%in%c(4,7))] <- "group47"

tf.dat.f <- tf.dat[which(rownames(tf.dat)%in%colnames(sub.dat)),]
all(rownames(tf.dat.f)==colnames(sub.dat))

# tf.dat.f$group <- sub.dat$group
# tmp.res <- aggregate(.~group,data=tf.dat.f,FUN=mean)
# rownames(tmp.res ) <- tmp.res$group
# tmp.res$group <- NULL

# caculate score and pvalue 
# 2021-6-1
#===================================================================================================
idx.09 <- which(sub.dat$group=="group09")
idx.16 <- which(sub.dat$group=="group16")
idx.47 <- which(sub.dat$group=="group47")

tmp.res <- t(apply(tf.dat.f,2,function(x){
    fc09 <- mean(x[idx.09])/mean(x[-idx.09])
    pvalue09 <- wilcox.test(x[idx.09],x[-idx.09])$p.value

    fc16 <- mean(x[idx.16])/mean(x[-idx.16])
    pvalue16 <- wilcox.test(x[idx.16],x[-idx.16])$p.value

    fc47 <- mean(x[idx.47])/mean(x[-idx.47])
    pvalue47 <- wilcox.test(x[idx.47],x[-idx.47])$p.value

    c(fc09,pvalue09,fc16,pvalue16,fc47,pvalue47)

}))
tmp.res <- data.frame(tmp.res)
colnames(tmp.res) <- c("FC09","pvalue09","FC16","pvalue16","FC47","pvalue47")



# caculate DEG 
#====================================================================================================
DefaultAssay(sub.dat) <- "RNA"
sub.dat@active.ident <- factor(sub.dat$group)
tmp.gene <- FindAllMarkers(sub.dat,assay="RNA",logfc.threshold=0.1)
tmp.gene <- tmp.gene[order(tmp.gene$avg_logFC,decreasing=T),]

deg <- tmp.gene[which(tmp.gene$avg_logFC>0&tmp.gene$p_val_adj<0.01),]

deg09 <- deg[which(deg$cluster=="group09"&deg$avg_logFC>1),]
write.table(as.vector(deg09$gene),file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Group09_DEG.txt",sep="\t",row.names=F,col.names=F,quote=F)

# the same cut off do not have DEGs 
deg16 <- deg[which(deg$cluster=="group16"&deg$avg_logFC>1),]
write.table(as.vector(deg16$gene),file="./tmp.txt",sep="\t",row.names=F,col.names=F,quote=F)

deg47 <- deg[which(deg$cluster=="group47"&deg$avg_logFC>1),]
write.table(as.vector(deg47$gene),file="./tmp.txt",sep="\t",row.names=F,col.names=F,quote=F)











###################################################################################################################################
# 2021-6-1
# plot result 
#==================================================================================================================================
# 0. TSNE plot show different group 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
dat$group <- "non"
dat$group[which(dat$seurat_clusters%in%c(0,9))] <- "group09"
dat$group[which(dat$seurat_clusters%in%c(1,6))] <- "group16"
dat$group[which(dat$seurat_clusters%in%c(4,7))] <- "group47"
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/sub.Tumor.DimPlot.pdf",useDingbats=F)
DimPlot(dat,group.by="group",cols=c("#2baf2b","#00a5dc","#ffc719","#ced7df"))
dev.off()




# 1. HPSC 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
# dat$group <- "group09"
# dat$group[which(dat$seurat_clusters%in%c(1,6))] <- "group16"
# dat$group[which(dat$seurat_clusters%in%c(4,7))] <- "group47"
# saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")

mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/sub.Tumor.HPSC.RDS")
all(colnames(dat)==rownames(mod))

mod$group <- factor(dat$group,levels=c("group47","group09","group16"))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/subtumor_HPSC.pdf",useDingbats=F)
boxplot(HPSC_C5_norm~group,data=mod,FUN=median,outline=F,ylim=c(0,1))
dev.off()






# 2.DEG analysis 
# find Markers to get DEG 
# then use enrichr to do hallmarks analysis 
#==============================================================================================================================
dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/MSigDB_Hallmark_2020_group09.txt",sep="\t",header=T)

dat$logOR <- log2(dat$Odds.Ratio)
dat <- dat[order(dat$Odds.Ratio,decreasing=T),]
dat <- dat[1:10,]
dat$Term <- factor(dat$Term,levels=rev(as.vector(dat$Term)))
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/sub.Tumor.Top10.hallmark.group09.pdf",useDingbats=F)
ggplot(dat,aes(x=logOR,y=Term)) + geom_point(aes(size= -log10(Adjusted.P.value),color= logOR))+
    scale_colour_gradientn(colours =  c("#598c14","#b5c327","#edd812","#faae40","#ff3c41","red"),values = c(0,0.4,0.5,0.8,1.0))+ 
    theme_bw()
dev.off()









# 3. heatmap show Immune gene  
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
dat@active.ident <- factor(dat$group)
tmp <- AverageExpression(dat,assay="RNA",feature=c("IDO1","CD274","KYNU","CD47","CD96","LGALS9"))$RNA
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/sub.Tumor.immune_heatmap.pdf",useDingbats=F)
pheatmap::pheatmap(tmp,scale="row",color=colorRampPalette(c('steelblue','white','red'))(50))
dev.off()







# 4. PySCENIC  show group16 expression Brain tissue TFs 
#===========================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
##################################################################################################################
# 2021-8-9
# use LCBM sub tumor RDS to calculate TFs
# tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM_subtumor/step3.auc_mtx.csv",sep=",",header=T) # do not OK
tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM/step3.auc_mtx.csv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
tf.dat.f <- tf.dat[which(rownames(tf.dat)%in%colnames(dat)),]

all(rownames(tf.dat.f)==rownames(dat@meta.data))

# analysis 
#===========================================================================================================================
idx.09 <- which(dat$group=="group09")
idx.16 <- which(dat$group=="group16")
idx.47 <- which(dat$group=="group47")

tmp.res <- t(apply(tf.dat.f,2,function(x){
    fc09 <- mean(x[idx.09])/mean(x[-idx.09])
    pvalue09 <- wilcox.test(x[idx.09],x[-idx.09])$p.value

    fc16 <- mean(x[idx.16])/mean(x[-idx.16])
    pvalue16 <- wilcox.test(x[idx.16],x[-idx.16])$p.value

    fc47 <- mean(x[idx.47])/mean(x[-idx.47])
    pvalue47 <- wilcox.test(x[idx.47],x[-idx.47])$p.value

    c(fc09,pvalue09,fc16,pvalue16,fc47,pvalue47)

}))
tmp.res <- data.frame(tmp.res)
colnames(tmp.res) <- c("FC09","pvalue09","FC16","pvalue16","FC47","pvalue47")

# caculate C1 and C6 TFs and plot result 
##################################################################################
tmp <- tmp.res[,c("FC16","pvalue16"),drop=F]
rs.f <- tmp[which(tmp$pvalue16<0.05),]
rs.f <- rs.f[order(rs.f$FC16,decreasing=T),]
plot.rs <- rs.f[1:10,]
plot.rs$gene <- rownames(plot.rs)
plot.rs$log10p <- -log10(plot.rs$pvalue16)
plot.rs$log10p[which(is.infinite(plot.rs$log10p))] <- max(plot.rs$log10p[-which(is.infinite(plot.rs$log10p))])
plot.rs$gene <- factor(plot.rs$gene,levels=rownames(plot.rs))
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/sub.Tumor.group16.TF.pdf",useDingbats=F)
ggplot(plot.rs,aes(x= gene ,y= FC16 ,fill= log10p)) + geom_bar(stat="identity",position="dodge") + theme_bw() +
    labs(x="Transcription factors",y="Fold change")

dev.off()






# tmp.res$group <- factor(tmp.res$group,levels=c("group09","group16","group47","Unknow"))




















#################################################################################
# USE https://www.proteinatlas.org/ENSG00000125820-NKX2-2/tissue
# show NKX2.2 is high expression in Brain 
#================================================================================




##################################################################################
# metabolism 
#=================================================================================

tmp <- data.table::fread("/public/workspace/lily/metastasis/data/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",)
tmp <- as.data.frame(tmp)
tmp.f <- tmp[,as.numeric(c(1,2,which(colnames(tmp)%in%as.vector(ann[which(ann$SMTS%in%c("Brain","Lung")),"SAMPID"]))))]
tmp.f$Name <- NULL
tmp.res <- aggregate(.~Description,data=tmp.f,FUN=median) # aggregate gene expression

rownames(tmp.res) <- tmp.res$Description
tmp.res$Description <- NULL

ann.f <- ann[which(ann$SAMPID%in%colnames(tmp.res)),]
rownames(ann.f) <- ann.f$SAMPID
saveRDS(tmp.res,file="/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung.RDS")
saveRDS(ann.f,file="/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung_ann.RDS")


#==================================================================================
library(GSEABase)
library(GSVA)
library(parallel)
genelist <- getGmt("/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt")
dat <- readRDS("/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung.RDS")
dat <- log2(dat+1)
result_gsva <- gsva(as.matrix(dat), genelist, kcdf="Gaussian", method="ssgsea",ssgsea.norm=TRUE, parallel.sz=10)

saveRDS(result_gsva,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GTEx_Lung_Brain_metabolism_GSVA.RDS")


###################################################################################
# 2021-6-3
# 2021-6-4 use ssGSEA to calculate metabolism mod 
#==================================================================================
# dat <- readRDS("/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GTEx_Lung_Brain_metabolism_ssGSEA.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung_ann.RDS")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
# mod.rs <- mod.analyze2(as.matrix(dat),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
# mod.rs <- as.data.frame(mod.rs)

# mod <- data.frame(mod.rs[,86:170])
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GTEx_Lung_Brain_metabolism_ssGSEA.RDS")


rownames(mod) <- gsub("\\.","-",rownames(mod))
all(rownames(mod)==rownames(ann))
mod$type <- ann$SMTS
tmp.res <- aggregate(.~type,data=mod,FUN=median)
rownames(tmp.res) <- tmp.res$type
tmp.res$type <- NULL
tmp.res <- data.frame(t(tmp.res))
tmp.res$diff <- tmp.res$Brain -tmp.res$Lung

tmp.res <- tmp.res[order(tmp.res$diff,decreasing=T),]

###############################################################################################################################################################
# use GTEX Brain and Lung data to show Brain specific metabolism pathway 
# Alanine__aspartate_and_glutamate_metabolism_norm
# Mannose_type_O.glycan_biosynthesis_norm
# Taurine_and_hypotaurine_metabolism_norm
# Butanoate_metabolism_norm
# Arginine_and_proline_metabolism_norm
# Oxidative_phosphorylation_norm
# and so on 


#============================================================================================================================
# 2021-6-4
# calculate metabolism for sub.Tumor 
#============================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
mod.rs <- mod.analyze2(as.matrix(dat[["RNA"]]@data),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)
saveRDS(mod.rs,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/sub.Tumor.metabolism.RDS")


#======================================================================================================================
# analysis result 
# 2021-6-5
#######################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/sub.Tumor.metabolism.RDS")
tmp.f <- data.frame(tmp.dat[,86:170])
tmp.f$group <- dat$group

tmp.res <- aggregate(.~group,data=tmp.f,FUN=median)
rownames(tmp.res) <- tmp.res$group
tmp.res$group <- NULL

tmp.res <- t(tmp.res)
res.f <- tmp.res[-which(rowSums(tmp.res)==0),]

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/sub.Tumor.metabolism.pdf",useDingbats=F,height=20,width=10)
pheatmap::pheatmap(res.f,scale="row")
dev.off()








#############################################################################################################################
# 2021-6-16
# calculate Nature Medicine signature 
#============================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/NatureMedicine/"))
mod <- mod.analyze2(as.matrix(dat[["RNA"]]@data),modlist,'/public/workspace/lily/MOD_file/NatureMedicine/',permN=0)
mod <- data.frame(mod)











# 2021-8-6
#===========================================================================================================
# run Pyscience to analysis subtype specific TFs 
# just use these cells
#===========================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
##################################################################################################################
# 2021-8-9
# use LCBM sub tumor RDS to calculate TFs
# tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM_subtumor/step3.auc_mtx.csv",sep=",",header=T) # do not OK
tf.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM/step3.auc_mtx.csv",sep=",",header=T)
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
tf.dat.f <- tf.dat[which(rownames(tf.dat)%in%colnames(dat)),]

all(rownames(tf.dat.f)==rownames(dat@meta.data))

# analysis 
#===========================================================================================================================
idx.09 <- which(dat$group=="group09")
idx.16 <- which(dat$group=="group16")
idx.47 <- which(dat$group=="group47")

tmp.res <- t(apply(tf.dat.f,2,function(x){
    fc09 <- mean(x[idx.09])/mean(x[-idx.09])
    pvalue09 <- wilcox.test(x[idx.09],x[-idx.09])$p.value

    fc16 <- mean(x[idx.16])/mean(x[-idx.16])
    pvalue16 <- wilcox.test(x[idx.16],x[-idx.16])$p.value

    fc47 <- mean(x[idx.47])/mean(x[-idx.47])
    pvalue47 <- wilcox.test(x[idx.47],x[-idx.47])$p.value

    c(fc09,pvalue09,fc16,pvalue16,fc47,pvalue47)

}))
tmp.res <- data.frame(tmp.res)
colnames(tmp.res) <- c("FC09","pvalue09","FC16","pvalue16","FC47","pvalue47")

#################################################################################################################
# 2021-8-9 
# plot 3D plot to show result group specific TFs 
#================================================================================================================
# another way to show TFs 
# add some info to plot result 
tmp.res$group <- "Unknow"
tmp.res$group[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC16,decreasing=T),],2)))] <- "group16"
tmp.res$group[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC47,decreasing=T),],2)))] <- "group47"
tmp.res$group[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC09,decreasing=T),],2)))] <- "group09"
tmp.res$group[which(rownames(tmp.res)=="CEBPD")] <- "specifc"
# set color 
colss <- forcats::fct_recode(as.factor(tmp.res$group),
                  "#2baf2b" = "group09","#00a5dc" = "group16",  # new <- old
                  "#ffc719" = "group47","#ced7df" = "Unknow","red" = "specifc"
    )


tmp.res$gene <- ""
tmp.res$gene[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC16,decreasing=T),],2)))] <- c("NKX2.2","SALL1")
# c("BHLHE41","NKX2.2","SALL1","ZNF384","ZNF589","ZNF770")
tmp.res$gene[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC47,decreasing=T),],2)))] <- c("E2F2","MYBL1")
# c("E2F2","E2F8","ETV4","MYBL1","NFYA","TEAD4")
tmp.res$gene[which(rownames(tmp.res) %in% rownames(head(tmp.res[order(tmp.res$FC09,decreasing=T),],2)))] <- c("NR4A2","ZNF880")
# c("ERG","NFIA","NR2F1","NR4A2","RUNX1","ZNF880")
tmp.res$gene[which(rownames(tmp.res)=="CEBPD")] <- "CEBPD"

library(plot3D)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/TF_specific_add.pdf",useDingbats=F)
with(tmp.res, scatter3D(x = FC09, y = FC47, z = FC16,
  pch = 21, cex = 1.5,col="black",bg=colss,
                   xlab = "Fold change group09",
                   ylab = "Fold change group47",
                   zlab = "Fold change group16", 
                   ticktype = "detailed",bty = "f",box = TRUE,
                   theta = -30, phi = 20, d=3,
                   colkey = FALSE))

with(tmp.res,text3D(x = FC09, y = FC47, z = FC16,col="black",bg=colss,
       colkey = FALSE, add = TRUE, 
       labels = tmp.res$gene))

dev.off()











###############################################################################################################################################################
# 2021-8-9
# try to use GSE110495 to analysis 
# 2021-9-15 use GSE110495 to plot results 
# check origin LT and BMIT 
#==============================================================================================================================================================
load("~/metastasis/data/verify/GSE110495/GSE110495.RData")
load("~/metastasis/data/verify/GSE110495/GSE110495_anno.RData")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(gse110495),c("Lung_gene","Brain_gene","BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$group <- anno[,2]
mod.f <- mod[which(mod$group%in%c("original","LT","BMIT")),]

# reshape 
library(reshape2)
library(ggplot2)
tmp.res <- reshape2::melt(mod.f[,c(4:6,10)])
colnames(tmp.res) <- c("group","type","score")
tmp.res$group <- factor(tmp.res$group,levels=c("original","LT","BMIT"))
tmp.res$type <- factor(tmp.res$type,levels=c("BMS_update_norm","Lung_gene_norm","Brain_gene_norm"))


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/GSE110495_3score_boxplot.pdf",useDingbats=F)
ggplot(tmp.res,aes(x=group,y=score,fill=type)) + geom_boxplot(outlier.color=NA) + scale_fill_manual(values=c("#ffc719","#2baf2b","#00a5dc"))
dev.off()















#==============================================================================================================================================================
# 2021-8-12
# need to make sure cell numbers of different cell type 
###############################################################################################################################################
# 0. self data 
libray(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
table(dat$celltype,dat$orig.ident)


# 1. D0927 and E927
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/D0927.RDS")
#check <- function(dat){

    DefaultAssay(dat) <- "RNA"
    pdf("~/tmp.pdf",width=10,height=10)
    DimPlot(dat,label=T)
    FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,order=F) # T cell
    FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','FCGR1A'),label=T,order=F) # Myeloid
    FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,order=T) # B cell 
    FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,order=T) # Oligodendrocyte
    FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,order=T) # Fibroblast/Vascular
    FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,order=T) # Endothelial
    FeaturePlot(dat,features=c("EGFR","EPCAM","PTPRC"),label=T,order=T,cols=c("lightgrey", "red"))

    VlnPlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),pt.size=0)
    VlnPlot(dat,features=c('CD19','CD68','FCGR3A','FCGR1A'),pt.size=0)
    VlnPlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),pt.size=0)
    VlnPlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),pt.size=0)
    VlnPlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),pt.size=0)
    VlnPlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),pt.size=0)
    VlnPlot(dat,features=c("EGFR","EPCAM","PTPRC"),pt.size=0)
    dev.off()

#}
# rm(list=ls())
# run again
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/E0927.RDS")


tmp <- matrix(c(58,37,136,83,192,256,488,1543,1021,90,770,421,101,32,1631,492,161,122,991,350,157,49,708,136),ncol=4,byrow=T)
colnames(tmp) <- c("A20190305","A20190312","T-Bsc1","E0927")
rownames(tmp) <- c("Endothelial","Fibroblast","group09","group16","group47","Oligodendrocyte")


1021	90	770	421
101	32	1631	492
161	122	991	350
157	49	708	136







































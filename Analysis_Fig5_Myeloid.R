
# 2022-6-30
# analysis about myeloid cells
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
Myeloid <- subset(dat,cells=which(dat$celltype.refine=="Myeloid"))

RE_INTE <- function(dat,sample){
    inte.list <- list()
    samplelist <- unique(dat@meta.data[,sample])
    for(i in 1:length(samplelist)){
        tmp <- subset(dat,cells=which(dat@meta.data[,sample]==samplelist[i]))
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

    return(inte)
}
Myeloid <- RE_INTE(Myeloid,sample="orig.ident")

saveRDS(Myeloid,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")






# Prepare for cell type classify 
#==========================================================================================================================
DefaultAssay(Myeloid) <- "RNA"
future::plan(multisession, workers=20)
marker <- FindAllMarkers(Myeloid,assay="RNA",logfc.threshold=0,min.pct=0.05)
library(SingleR)
ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_NovershternHematopoieticData.RDS")

ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster1<-SingleR(test=as.matrix(Myeloid@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=Myeloid$seurat_clusters,method="cluster")


Myeloid$celltype.refine <- "unclassify"
Myeloid$celltype.refine[which(Myeloid$seurat_clusters%in%c(0,1,3,5,6,8,11,12,13,14,15,20))] <- "Macrophage"
Myeloid$celltype.refine[which(Myeloid$seurat_clusters%in%c(4))] <- "Alveolar Mac"
Myeloid$celltype.refine[which(Myeloid$seurat_clusters%in%c(2,7,9,23,17))] <- "Monocyte"
Myeloid$celltype.refine[which(Myeloid$seurat_clusters%in%c(10))] <- "Mast"
Myeloid$celltype.refine[which(Myeloid$seurat_clusters%in%c(16,21,22))] <- "DC"

saveRDS(Myeloid,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")





# 2022-8-22
# changed cell type and re-plot result 
# 0. plot tsne-result 
#============================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/DimPlot_Myeloid.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype.refine",cols=c("#7552cc","#ea4c89","#a25016","#00334e","#0077c8","#aea400","#bfbfbf"),reduction="tsne",raster=F)
dev.off()




# 1. plot percentage 
tmp <- table(dat$type_group,dat$celltype.refine)
tmp.f <- tmp[,-grep("unclassify",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("nMLUAD","MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Alveolar Mac","DC","Macrophage","Mast","MG","Monocyte"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/Myeloid_subtype_percent.pdf",useDingbats=F)
cols <- c("#7552cc","#ea4c89","#F29403","#00334e","#0077c8","#aea400")
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



# # 2. pro-inflammatory and anti-inflammatroy factor 
# gene <- c("IL1B", "IL6", "IFNG","IL12A","NOS2",   "IL4", "IL13", "ARG1", "CXCL8","IDO1","IL10","TGFB1","VEGFA")

# percent_feature <- function(dat,genelist,group){
#     res.list <- c()
#     for(i in 1:length(genelist)){
#         dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
#         if(all(dat$tmp_gene=="N")){
#            res.list[[i]] <- rep(0,length=length(table(dat@meta.data[,group])))
#         }else{
#            res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
#         }
        
#         names(res.list)[i] <- genelist[i]
#     }
#     return(res.list)
# }
# dat$celltype_group <- paste0(dat$celltype.refine,"_",dat$type_group)
# res.dat.list <- percent_feature(dat,gene,group="celltype_group")

# list2mat <- function(res.list){
#     res.dat <- matrix(unlist(res.list),ncol=length(res.list[[1]]),byrow=T)
#     rownames(res.dat) <- names(res.list)
#     colnames(res.dat) <- names(res.list[[1]])
#     return(res.dat)
# }

# mat.res <- list2mat(res.dat.list)

# pheatmap::pheatmap(mat.res[,7:9],cluster_rows=F,cluster_cols=F,scale="row")



# # find DEG
# mac <- subset(dat,cells=which(dat$celltype.refine=="Macrophage"))
# DefaultAssay(mac) <- "RNA"
# mac@active.ident <- factor(mac$type_group)
# macmarker <- FindAllMarkers(mac,logfc.threshold=0,min.pct=0.05)

# tmp.res <- macmarker[which(macmarker$cluster=="LCBM"&macmarker$p_val_adj<0.05),]
# tmp.res$class <- "unknow"
# tmp.res$class[which(tmp.res$avg_logFC>0.5 &tmp.res$pct.1>0.05)] <- "up"
# tmp.res$class[which(tmp.res$avg_logFC< -0.5 &tmp.res$pct.1>0.05)] <- "dn"
# tmp.res$class <- factor(tmp.res$class,levels=c("dn","unknow","up"))

# # add some label
# tmp.res$label <- NA
# tmp.res$label[which(tmp.res$gene=="CCL20")] <- "CCL20"
# tmp.res$label[which(tmp.res$gene=="FN1")] <- "FN1"
# tmp.res$label[which(tmp.res$gene=="MMP12")] <- "MMP12"
# tmp.res$label[which(tmp.res$gene=="CD163")] <- "CD163"




# library(ggplot2)
# library(ggrepel)

# pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/Myeloid_LCBM_DEG.pdf",useDingbats=F)
# ggplot(tmp.res,aes(x=avg_logFC,y=pct.1,group=class,color=class)) + geom_point(size=2.5) + 
#     scale_colour_manual(values=c("#2e9df7","#d7d7d8","#ec2c22")) +  geom_vline(aes(xintercept= -0.5),color="#2e9df7") +
#     geom_vline(aes(xintercept= 0.5),color="#ec2c22") + geom_hline(aes(yintercept= 0.05),color="black")+
#     geom_text_repel(aes(x=avg_logFC,y=pct.1,label=label))

# dev.off()





#===========================================================================================================================
# calculate macrophage MDM and MG
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/inte_16S_Myeloid_MDM_MG_mod.RDS")


#==================================================================================================================
# 2022-8-5 
# use marker gene to identify MG and MDM
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
FeaturePlot(dat,features=c("SLC2A5","P2RY12","NAV3"),label=T,label.size=6)
dat$celltype.refine[which(dat$seurat_clusters==8)] <- "MG"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")



#=================================================================================================================================================
# check OSM expression percentage 
dat$OSM <- ifelse(dat[["RNA"]]@data["OSM",]>0,"Y","N")

library(ggpubr)
#OSM
tmp.dat <- apply(table(dat$OSM,dat$celltype.refine),1,function(x){x/sum(x)})
tmp.dat <- data.frame(tmp.dat)
tmp.dat$celltype <- rownames(tmp.dat)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/Myeloid_OSM_celltype_pie.pdf",useDingbats=F)
ggpie(tmp.dat, 'Y',  
    fill = 'celltype',
    palette = c("#7552cc","#ea4c89","#F29403","#00334e","#0077c8","#aea400","#bfbfbf"), 
    label = paste0(round(tmp.dat$Y * 100,1),"%"), 
    lab.pos = 'in', lab.font = c(4, 'white') #设置标签，标签的位置在图的内部，标签的大小为4， 颜色为白色.
) 
dev.off()




# plot OSM+ macrophage in type groups 
subdat <- subset(dat,cells=which(dat$celltype.refine=="Macrophage"&dat$OSM=="Y"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/OSM+MAC_type_group.pdf",useDingbats=F)
barplot(table(subdat$type_group)/ncol(subdat),ylim=c(0,0.5))
dev.off()






#================================================================================================================================================
# 2022-8-22
# 接下来就应该聚焦Macrophage
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Macrophage"))
DefaultAssay(dat) <- "RNA"


# find DEG
dat@active.ident <- factor(dat$type_group)
future::plan("multisession", workers=20)
options(future.globals.maxSize=800000000)
macmarker <- FindAllMarkers(dat,logfc.threshold=0,min.pct=0.05)




# 1. Diffusion map
#=================================================================================================================================================
library(Seurat)
library(destiny)
library(SingleCellExperiment)
# Diffusion map for All type group
# tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
# dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine%in%c("Macrophage","Monocyte")))
# sce <- as.SingleCellExperiment(dat)
# dm <- DiffusionMap(sce, verbose = TRUE)
# saveRDS(dm,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Monocyte_Macrophage_Diffusionmap.RDS")


# Diffusion map for LCBM samples
# tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
# dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine%in%c("Macrophage","Monocyte")& tmp.dat$type_group=="LCBM"))
# sce <- as.SingleCellExperiment(dat)
# dm <- DiffusionMap(sce, verbose = TRUE)
# saveRDS(dm,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/LCBM_Monocyte_Macrophage_Diffusionmap.RDS")

tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine%in%c("Macrophage","Monocyte")& tmp.dat$type_group=="LCBM"))
dm <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/LCBM_Monocyte_Macrophage_Diffusionmap.RDS")
dat$OSM <- ifelse(dat[["RNA"]]@data["OSM",]>0,"Y","N")
dat$OSM_value <- as.numeric(dat[["RNA"]]@data["OSM",])


library(ggplot2)
meta<-dat@meta.data
meta$id<-rownames(meta)
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1],
                  DC2 = eigenvectors(dm)[, 2],
                  DC3 = eigenvectors(dm)[, 3],
                  DC4 = eigenvectors(dm)[, 4],
      id = names(optimal_sigma(dm)))
dest<-merge(meta,tmp)      


# dest$order <- ifelse(dest$OSM=="P",1,2)

library(dplyr)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/OSM_Diffusionmap.pdf",useDingbats=F)
ggplot(dest %>% arrange(OSM), aes(x = DC1, y = DC2, colour = OSM)) +
  geom_point()+ 
  scale_color_manual(values=c('#e9e8dd','#9fbb58'))+
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
dev.off()


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/Celltype_Diffusionmap.pdf",useDingbats=F)
ggplot(dest, aes(x = DC1, y = DC2, colour = celltype.refine)) +
  geom_point()+ 
  scale_color_manual(values=c("#F29403","#aea400"))+
  xlab("Diffusion component 1") + 
  ylab("Diffusion component 2") +
  theme_classic()
dev.off()



# 2. Find DEG 
# DEG in OSM+ vs. OSM-
#==========================================================================================================================
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Macrophage"))
DefaultAssay(dat) <- "RNA"
subdat <- subset(dat,cells=which(dat$type_group=="LCBM"))

subdat$OSM <- ifelse(subdat[["RNA"]]@data["OSM",]>0,"P","N")
subdat@active.ident <- factor(subdat$OSM)
macmarker <- FindAllMarkers(subdat,logfc.threshold=0,min.pct=0.05)

saveRDS(macmarker,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/LCBM_MacOSM_DEG.RDS")
macmarker <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/LCBM_MacOSM_DEG.RDS")
macmarker$pct.diff <- abs(macmarker$pct.1 - macmarker$pct.2)
plotdat <- macmarker[which(macmarker$cluster=="P"&macmarker$p_val_adj<0.05),]

plotdat$group <- "unknow"
plotdat$group[which(plotdat$avg_logFC >0.2&plotdat$pct.diff>0.2)] <- "up"

# add gene name

plotdat$geneName <- plotdat$gene
plotdat$geneName[-which(plotdat$geneName%in%c("OSM","HBEGF","NR4A3","NR4A2","NR4A1"))] <- NA


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/OSM+Macrophage_DEG.pdf",useDingbats=F)

library(ggplot2)
library(ggrepel)
ggplot(plotdat[which(plotdat$avg_logFC>0),],aes(x=pct.diff,y=avg_logFC,group=group,color=group)) + 
    geom_point(size=3) + 
    scale_colour_manual(values=c("#d7d7d8","#ec2c22")) +
    geom_label(aes(label = geneName), size = 3)+ labs(title="Up-regulation in OSM+ Macrophage") +
    geom_vline(aes(xintercept= 0.2),color="red") +
    geom_hline(aes(yintercept= 0.2),color="red") 

dev.off()





# 3. HBEGF in LCBM 
# GSE200563 Rho=0.39, P= 0.06
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM")),]
dat <- tmp.dat[,rownames(sampleinfo)]


# GSE14108
# GSE14108 data 
# Rho = 0.75 p < 0.001
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
cor.test(as.numeric(dat["OSM",]),as.numeric(dat["HBEGF",]),method="spearman")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/GSE14108_OSM_HBEGF_cor.pdf",useDingbats=F)
plot(as.numeric(dat["OSM",]),as.numeric(dat["HBEGF",]),main=" GSE14108 ")
abline(lm(as.numeric(dat["HBEGF",])~as.numeric(dat["OSM",])),col="red")
legend("topright",legend=paste0("rho=",cor.test(as.numeric(dat["OSM",]),as.numeric(dat["HBEGF",]),method="spearman")$estimate,
    " P =",round(cor.test(as.numeric(dat["OSM",]),as.numeric(dat["HBEGF",]),method="spearman")$p.value,2))
)

dev.off()




















# 4. HBEGF in small large Brain metastasis
# GSE137762
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
Hbegf.res <- as.data.frame(dat["Hbegf",])
Hbegf.res$type <- "unknow"
colnames(Hbegf.res)[1] <- "exp"
Hbegf.res$type <- gsub("[0-9]$","",rownames(Hbegf.res))

# Osm
Osm.res <- as.data.frame(dat["Osm",])
Osm.res$type <- "unknow"
colnames(Osm.res)[1] <- "exp"
Osm.res$type <- gsub("[0-9]$","",rownames(Osm.res))
tmp.f <- Osm.res[grep("BMDM|Blood",rownames(Osm.res)),]
tmp.f <- tmp.f[1:14,]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM"))

data <- data.frame(aggregate(exp~type,data=tmp.f,FUN=mean))
data$sd <- aggregate(exp~type,data=tmp.f,FUN=sd)[,2]
data$se = data$sd / sqrt(5)  #sd(vec) / sqrt(length(vec))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/GSE137762_Osm_barplot.pdf",useDingbats=F)
ggplot(data) +
    geom_bar( aes(x=type, y=exp), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=type, ymin=exp-se, ymax=exp+se), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme_bw()
dev.off()


#=======================================================================================================================
# just show Hbegf not show 
#=======================================================================================================================
tmp.f <- Hbegf.res[grep("BMDM|Blood",rownames(Hbegf.res)),]
tmp.f <- tmp.f[1:14,]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM"))

data <- data.frame(aggregate(exp~type,data=tmp.f,FUN=mean))
data$sd <- aggregate(exp~type,data=tmp.f,FUN=sd)[,2]
data$se = data$sd / sqrt(5)  #sd(vec) / sqrt(length(vec))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/GSE137762_Hbegf_barplot.pdf",useDingbats=F)
ggplot(data) +
    geom_bar( aes(x=type, y=exp), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=type, ymin=exp-se, ymax=exp+se), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme_bw()
dev.off()










# 5. NR4A familiy with NFKB
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)

cor.test(as.numeric(mod[,2]),as.numeric(dat["NR4A2",]),method="spearman")


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/GSE14108_NR4A2_NFKB_cor.pdf",useDingbats=F)

plot(as.numeric(mod[,2]),as.numeric(dat["NR4A2",]),main=" GSE14108 ")
abline(lm(as.numeric(dat["NR4A2",])~as.numeric(mod[,2])),col="red")
legend("topright",legend=paste0("rho=",cor.test(as.numeric(mod[,2]),as.numeric(dat["NR4A2",]),method="spearman")$estimate,
    " P =",round(cor.test(as.numeric(mod[,2]),as.numeric(dat["NR4A2",]),method="spearman")$p.value,2))
)

dev.off()





# 6. Immune suppression gene 
# https://www.cell.com/cell-reports/pdf/S2211-1247(20)31560-6.pdf


percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
           res.list[[i]] <- rep(0,length=length(table(dat@meta.data[,group])))
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}


tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
dat <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Macrophage"))
DefaultAssay(dat) <- "RNA"
subdat <- subset(dat,cells=which(dat$type_group=="LCBM"))

subdat$OSM <- ifelse(subdat[["RNA"]]@data["OSM",]>0,"P","N")
subdat@active.ident <- factor(subdat$OSM)

features <- c("ARG1","PTGS2","KYNU","QPRT",
"IDO1","LGALS9",
"PVR","CD274","CD80","BTLA"
)

tmp.list <- percent_feature(subdat,genelist=features,group="OSM")

data <- data.frame(row.names=names(tmp.list),gene=names(tmp.list))
data$pctN <- as.numeric(sapply(tmp.list,function(x){as.numeric(scale(x,center=F))[1]}))
data$pctP <- as.numeric(sapply(tmp.list,function(x){as.numeric(scale(x,center=F))[2]}))

plotdat <- reshape2::melt(data,var="gene")
colnames(plotdat) <- c("gene","group","percent")
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_immunesuppresion.pdf",useDingbats=F)
ggplot(data=plotdat,aes(x=group,y=gene))+geom_point(aes(color=percent,size=I(7)))+
    scale_color_gradientn(colors=c("#b9cfed","#f4f7f9","#ee4f4f"))+
    theme_bw()
dev.off()








# 2022-9-12
# GSEA analysis gene 
macmarker <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/LCBM_MacOSM_DEG.RDS")
plotdat <- macmarker[which(macmarker$cluster=="P"&macmarker$p_val_adj<0.05),]

# get gene set
# https://www.frontiersin.org/articles/10.3389/fimmu.2019.01462/full
# M1 : CD80, CD86, CIITA,MHC-II,COX-2, iNOS,TNF-α, IL1-β, IL-6, IL-12, and IL-23,STAT3,STAT1,IRF4,AP1
# M2 : CD206, IL-1R, TGF-β, IGF-1,ARG1,IL-8,VEGF-A,CD163, CD36,STAT6, GATA3,PPARγ,FIZZ1,PDGF

M1geneset <- c("CD80","CD86","CIITA","HLA-DRA1","HLA-DRB1","HLA-DRB3","HLA-DPB1",
    "HLA-DPA1","HLA-DQB1","COX2","NOS2","TNFA","IL1B","IL6","IL12A","IL12B","IL23A",
    "STAT3","STAT1","IRF4","FOS","JUN")
M2geneset <- c("CD206","IL1R1","TGFB1","IGF1","ARG1","IL8","VEGFA","CD163","CD36","STAT6","GATA3","PPARG","FIZ1","PDGFA")
act <- c("TNFSF4","CD80","TNFSF15","ICOSLG","CD40","CD86","TNFSF18")
inb <- c("PDCD1LG2","LGALS9","BTNL2","HHLA2","CD276","VTCN1","CD274","IDO1","IDO2")

# analysis 
library(clusterProfiler)
library(enrichplot)

plotdat <- plotdat[order(plotdat$avg_logFC,decreasing=T),]
M1.df  <-  as.data.frame(M1geneset)
colnames(M1.df) <- NULL
M1.df$term <- "M1"
M2.df <- as.data.frame(M2geneset)
colnames(M2.df) <- NULL
M2.df$term <- "M2"
geneset <- rbind(M1.df,M2.df)
colnames(geneset) <- c("gene","term")
geneset <- geneset[,c("term","gene")] # col1 must be term 

genelist <- plotdat$avg_logFC
names(genelist) <- plotdat$gene
KEGG <- GSEA(geneList=genelist,TERM2GENE = geneset) #GSEA分析


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_M1_M2_GSEA.pdf",useDingbats=F)
gseaplot2(KEGG, "M1", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T)
gseaplot2(KEGG, "M2", color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T)
dev.off()


#==================================================================================================================================
# for gmt files 
# analysis 

plotdat <- plotdat[order(plotdat$avg_logFC,decreasing=T),]
genelist <- plotdat$avg_logFC
names(genelist) <- plotdat$gene

# NES 2.17
gmt <- read.gmt("/public/workspace/lily/MOD_file/GOBP_INFLAMMATORY_RESPONSE.v2022.1.Hs.gmt")
KEGG <- GSEA(geneList=genelist,TERM2GENE = gmt,pvalueCutoff=2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_inflammtory_GSEA.pdf",useDingbats=F)
gseaplot2(KEGG, 1, color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T)
dev.off()


# NES -1.39
gmt <- read.gmt("/public/workspace/lily/MOD_file/GOBP_PHAGOCYTOSIS.v2022.1.Hs.gmt")
KEGG <- GSEA(geneList=genelist,TERM2GENE = gmt,pvalueCutoff=2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_phagocytosis_GSEA.pdf",useDingbats=F)
gseaplot2(KEGG, 1, color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T)
dev.off()

# NES 1.47
gmt <- read.gmt("/public/workspace/lily/MOD_file/BIOCARTA_NFKB_PATHWAY.v2022.1.Hs.gmt")
KEGG <- GSEA(geneList=genelist,TERM2GENE = gmt,pvalueCutoff=2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_NFKB_GSEA.pdf",useDingbats=F)
gseaplot2(KEGG, 1, color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T)
dev.off()


# NES -1.47
gmt <- read.gmt("/public/workspace/lily/MOD_file/GOBP_ACTIVATION_OF_IMMUNE_RESPONSE.v2022.1.Hs.gmt")
KEGG <- GSEA(geneList=genelist,TERM2GENE = gmt,pvalueCutoff=2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/scRNA_LCBM_Mac_OSM_Immuneact_GSEA.pdf",useDingbats=F)
gseaplot2(KEGG, 1, color = "firebrick", rel_heights=c(1, .2, .6),pvalue_table=T) 
dev.off()









# 2022-9-15
# make a immun suppression mod 
# # https://www.cell.com/cell-reports/pdf/S2211-1247(20)31560-6.pdf
gene <- c("ARG1","PTGS2","KYNU","QPRT",
"IDO1","LGALS9",
"PVR","CD274","CD80","BTLA"
)

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"immunosuppression",out=paste0("/public/workspace/lily/Lung2Brain/Version6/","immunosuppression",".mod")) # make a mod file 

# Rho = 0.308 p =0.11
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")






tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM")),]
dat <- tmp.dat[,rownames(sampleinfo)]

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("immunosuppression"),"/public/workspace/lily/Lung2Brain/Version6/",permN=0)
mod <- as.data.frame(mod)

cor.test(as.numeric(dat["OSM",]),as.numeric(mod[,2]),method="spearman")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig5/GSE200563_OSM_Immunosuppression_cor.pdf",useDingbats=F)
plot(as.numeric(mod[,2]),as.numeric(dat["OSM",]),main=" GSE200563 ")
abline(lm(as.numeric(dat["OSM",])~as.numeric(mod[,2])),col="red")
legend("topright",legend=paste0("rho=",cor.test(as.numeric(mod[,2]),as.numeric(dat["OSM",]),method="spearman")$estimate,
    " P =",round(cor.test(as.numeric(mod[,2]),as.numeric(dat["OSM",]),method="spearman")$p.value,2))
)
dev.off()









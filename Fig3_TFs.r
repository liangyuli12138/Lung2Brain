
# 2021-11-17 
# this program is used to get CEBPD and CCL20 
# 0. BMS hallmark analysis 
# 1. maoding genes hallmark analysis 
# 2. 4 TFs activated genes and gene enrichment analysis
# 3. how about 脂代谢 and immune
# maybe need scRNA seq data ? 
# 4. cellchat for BMS high and BMS low interaction with TME 
#========================================================================================================================================

# 0. get maoding and plot 
# map gene 
dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
gene <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Signature/BMS_update_gene.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

tmp <- sort(apply(dat,1,function(x){cor.test(as.numeric(x),mod$BMS_update_norm,method="spearman")$estimate}))
tmp.f <- tmp[-which(names(tmp)%in%gene)] # do not get signature gene

write.table(names(tail(tmp.f,20)),file="~/tmp.txt",quote=F,col.names=F,row.names=F)

# for another data 
# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")

# now plot result
library(pheatmap)
mod <- mod[order(mod$BMS_update_norm,decreasing=T),]
dat.f <- dat[which(rownames(dat)%in%names(tail(tmp.f,20))),rownames(mod)]
dat.f <- dat.f[rev(names(tail(tmp.f,20))),]

# add row and column annotation
ann_row <- data.frame(row.names=rev(names(tail(tmp.f,20))),correlation=rev(unname(tail(tmp.f,20))))
ann_col <- data.frame(row.names=rownames(mod),BMS=mod$BMS_update_norm)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/GSE14108_maoding20gene.pdf",useDingbats=F)
pheatmap(dat.f,cluster_cols=F,cluster_rows=F,scale="row",annotation_row=ann_row ,
    annotation_col=ann_col,
    color = colorRampPalette(c("steelblue","#B9D0E2", "white", "#FF5252","red"))(100))
dev.off()




# 2021-11-22
# use BMS signature 
dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Data/EnrichR/MSigDB_Hallmark_2020_BMS.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log2(dat$Adjusted.P.value)
sub.dat <- dat[which(dat$Adjusted.P.value<0.05),]
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/BMS_hallmark.pdf",useDingbats=F)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+
    geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip() +
    scale_fill_gradientn(colors=c("lightgrey","steelblue"),limits=c(2, ceiling(max(sub.dat$logP))))
dev.off()








#============================================================================================================================================
# 1. plot GSE14108 hallmark result 
# Gene enrichment analysis 

dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Data/EnrichR/MSigDB_Hallmark_2020_GSE14108_20maodinggene.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log2(dat$Adjusted.P.value)
sub.dat <- head(dat,10)
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/GSE14108_maoding20gene_hallmark.pdf",useDingbats=F)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+
    geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip() +
    scale_fill_gradientn(colors=c("lightgrey","steelblue"))
dev.off()








#=============================================================================================================================================
# 2. plot sangji plot 
#=============================================================================================================================================
require(igraph)
require(rCharts)
require(rjson)
require(plyr)
require(reshape2)
require(networkD3)
options(stringsAsFactors=F)
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)
#=============================================================================================================================================
gene <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Signature/BMS_update_gene.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
tmp.res <- trrust.dat[which(trrust.dat$source%in%gene&trrust.dat$effect=="Activation"),]

# links 
links <- tmp.res[,c("source","traget","effect")]
links$effect <- 1
tmp <- to_lodes_form(links,axes = 1:2,id = "Cohort")


# plot 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/Sangjiplot_TF.pdf",useDingbats=F,height=10)
ggplot(tmp,aes(x =factor(x,level = c("source","traget")),y=effect,stratum = stratum, alluvium = Cohort,fill = stratum, label =stratum)) +
geom_flow( width = 1/4) + # flow
geom_stratum( width = 1/4,linetype=0,size=0.5,alpha =0.5,color = "black")+ 
geom_text(stat ="stratum" , size =3) + #添加名字
scale_x_discrete(limits = c() )+ #去掉横坐标轴
theme_bw() +
theme(legend.position="none") +
scale_fill_manual(values = colorSpace)
dev.off()








#=============================================================================================================================================
# 3. GSEA plot for gene set 
# 
#=============================================================================================================================================
# for CEBPD

dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Data/EnrichR/MSigDB_Hallmark_2020_CEBPD_act_gene.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log2(dat$Adjusted.P.value)
sub.dat <- dat[which(dat$Adjusted.P.value<0.05),]
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/CEBPD_act_.pdf",useDingbats=F)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+
    geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip() +
    scale_fill_gradientn(colors=c("lightgrey","steelblue","steelblue"),limits=c(2, ceiling(max(sub.dat$logP))))
dev.off()



# and for KLF5
dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Data/EnrichR/MSigDB_Hallmark_2020_KLF5_act_gene.txt",sep="\t",header=T)
library(ggplot2)
dat$logP <- -log2(dat$Adjusted.P.value)
sub.dat <- dat[which(dat$Adjusted.P.value<0.05),]
sub.dat <- sub.dat[order(sub.dat$logP,decreasing=T),]
sub.dat$Term <- factor(sub.dat$Term,levels=rev(sub.dat$Term))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/KLF5_act_.pdf",useDingbats=F)
ggplot(sub.dat,aes(x=Term,y=Odds.Ratio,fill=logP))+
    geom_bar(stat="identity",position = position_dodge(0.5),width=0.5) + coord_flip() +
    scale_fill_gradientn(colors=c("lightgrey","#F2BED7","#F2BED7"),limits=c(2, ceiling(max(sub.dat$logP))))
dev.off()










#=============================================================================================================================================
# 4. check LCBM have higher expression of CEBPD and KLF5 ?
# 
#=============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
dat@active.ident <- dat$seurat_clusters


dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")
sub.dat <- subset(dat,cells=which(dat$seurat_clusters==6))
sub.dat$group[which(sub.dat$group=="PRIMARY"&sub.dat$sample=="LX675")] <- "ADVANCED"

sub.dat@active.ident <- factor(sub.dat$group)





#============================================================================================================================================
# Caculate for scRNA-seq data 
# 0.25 percentage and 0.75 percentage for BMS group
# 2021-11-18
#============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

datname <- "pur_LCBM"
filepath <- "/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/"
# 0. define celltype
dat$BMS.group <- "inte"
dat$BMS.group[which(dat$BMS_update<quantile(dat$BMS_update,0.25))] <- "low"
dat$BMS.group[which(dat$BMS_update>quantile(dat$BMS_update,0.75))] <- "high"
DefaultAssay(dat) <- "RNA"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

# 1. plot CEBPD expression percentage 
pdf(paste0(filepath,datname,"_CEBPD_KLF5_featureplot.pdf"),useDingbats=F,width=18)
FeaturePlot(dat,features=c("CEBPD"),split.by="BMS.group")
FeaturePlot(dat,features=c("KLF5"),split.by="BMS.group")
dev.off()



# 2. plot metabolism 
# first run metabolism in APP
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/Metabolism/final_metabolism_pathways_activity_pur_LCBM_BMS.group.txt",header=T,sep="\t")
dat <- tmp.dat[c(59,75,1,3,16,69),c(2,3)]
dat$pathway <- rownames(dat)
library(ggplot2)
p.dat <- reshape2::melt(dat)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/pur_LCBM_metabolism.pdf",useDingbats=F)
ggplot(p.dat,aes(x=variable,y=pathway))+geom_point(aes(size=value,color=variable))
dev.off()

















#==========================================================================================================================================
# bulik data verify 
# CEBPD with apido,purity,and BMS 
# GSE14108 and E-MTAB is ok
#==========================================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GSE14108_estimate_score.gct",skip = 2,header = T)

# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
# scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/E_MTAB_estimate_score.gct",skip = 2,header = T)

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$purity <- as.numeric(as.vector(scores[4,-c(1,2)]))
mod$immune <- as.numeric(as.vector(scores[2,-c(1,2)]))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.rs <- mod.analyze2(as.matrix(dat),c("HALLMARK_ADIPOGENESIS","HALLMARK_FATTY_ACID_METABOLISM"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod.rs <- as.data.frame(mod.rs)

cor.test(as.numeric(dat["CEBPD",]),mod$purity,method="spearman")
cor.test(as.numeric(dat["CEBPD",]),mod$BMS_update_norm,method="spearman")
cor.test(as.numeric(dat["CEBPD",]),mod.rs$HALLMARK_ADIPOGENESIS_norm,method="spearman")
cor.test(as.numeric(dat["CEBPD",]),mod.rs$HALLMARK_FATTY_ACID_METABOLISM_norm,method="spearman")

# cor.test(mod$BMS_update_norm,mod.rs$HALLMARK_ADIPOGENESIS_norm,method="spearman")
# cor.test(mod$BMS_update_norm,mod.rs$HALLMARK_FATTY_ACID_METABOLISM_norm,method="spearman")

############################################################################################################################################
# plot result 

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/GSE14108_CEBPD_BMS_Fatty.pdf",useDingbats=F)
plot(as.numeric(dat["CEBPD",]),mod$purity,main="CEBPD with purity",xlab="CEBPD expression",ylab="estimate purity score",pch=19)
abline(lm(mod$purity~as.numeric(dat["CEBPD",])),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",-0.41," pvalue=",0.031))

plot(as.numeric(dat["CEBPD",]),mod$BMS_update_norm,main="CEBPD with BMS",xlab="CEBPD expression",ylab="BMS score",pch=19)
abline(lm(mod$BMS_update_norm~as.numeric(dat["CEBPD",])),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.34," pvalue=",0.08))

plot(as.numeric(dat["CEBPD",]),mod.rs$HALLMARK_ADIPOGENESIS_norm,main="CEBPD expression",xlab="HALLMARK_ADIPOGENESIS",ylab="estimate immune score",pch=19)
abline(lm(mod.rs$HALLMARK_ADIPOGENESIS_norm~as.numeric(dat["CEBPD",])),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.67," pvalue<",0.001))

plot(as.numeric(dat["CEBPD",]),mod.rs$HALLMARK_FATTY_ACID_METABOLISM_norm,main="CEBPD expression",xlab="HALLMARK_FATTY_ACID_METABOLISM",ylab="estimate immune score",pch=19)
abline(lm(mod.rs$HALLMARK_FATTY_ACID_METABOLISM_norm~as.numeric(dat["CEBPD",])),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.59," pvalue=",0.001))

dev.off()














#==========================================================================================================================================
# 4. run Cellchat for BMS high and BMS low 
library(Seurat)
library(CellChat)

LCBM.tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
LCBM.tumor <- subset(LCBM.tumor,cells=which(LCBM.tumor$BMS.group%in%c("high","low")))
ntumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/inte9_ntumor.RDS")
LCBM.ntumor <- subset(ntumor,cells=which(ntumor$type_group=="LCBM"))
LCBM.ntumor <- subset(LCBM.ntumor,cells=-which(LCBM.ntumor$type=="unknow"))
LCBM.ntumor$BMS.group <- LCBM.ntumor$type
dat <- merge(LCBM.tumor,LCBM.ntumor)

# Run CellChat 
data.input <- GetAssayData(dat,assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat@meta.data) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "BMS.group")
cellchat <- setIdent(cellchat, ident.use = "BMS.group") # set "labels" as default cell identity
levels(cellchat@idents)

CellChatDB <- CellChatDB.human
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel

#####################################################################################################
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (should use)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat) # 1W cells 5min
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

# saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.subtumor.CellChat.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/cellchat_LCBM_interaction.pdf",useDingbats=F)
groupSize <- table(cellchat@idents)
netVisual_circle(cellchat@net$count, sources.use = c(4,5),targets.use=c(1:3,6:8),vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$count, sources.use = 4,targets.use=c(1:3,6:8),vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$count, sources.use = 5,targets.use=c(1:3,6:8),vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()



pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3_TFs/cellchat_LCBM_interaction_num.pdf",useDingbats=F)
barplot(c(1915,440),names=c("BMS.high","BMS.low"),main="BMS.interaction.num",ylim=c(0,2000))
dev.off()














































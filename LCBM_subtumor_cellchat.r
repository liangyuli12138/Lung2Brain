
# 2021-7-12
# this program is used to analysis LCBM  different group cell cell communication with different cell type 
#============================================================================================================================================
##############################################################################################################################################################
# 2021-7-10 
# check if Oligodendrocyte expression 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
myeloid <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
tcell <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/inte9_ntumor.RDS")


myeloid.sub <- subset(myeloid,cells=which(myeloid$type_group=="LCBM"))
myeloid.sub$type <- "Myeloid"
tcell.sub <- subset(tcell,cells=which(tcell$type_group=="LCBM"))
tcell.sub$type <- "Tcell"
tmp.sub <- subset(tmp.dat,cells=which(tmp.dat$type%in%c("Endothelial","Fibroblast","Oligodendrocyte")&tmp.dat$type_group=="LCBM"))
dat$celltype <- dat$group
dat$type <- dat$group


tmp.sub$hbmarker <- NULL
tmp.sub$llymarker <- NULL
tmp.sub$celltype <- NULL
tmp.sub$putativeTumor3 <- NULL
tmp.sub$res.2 <- NULL
tmp.sub$Phase <- NULL
tmp.sub$G2M.Score <- NULL
tmp.sub$S.Score <- NULL
tmp.sub$percent.mito <- NULL
tmp.sub$nUMI <- NULL
tmp.sub$nGene <- NULL
tmp.sub$RNA_snn_res.2 <- NULL

# ready to 
tmp.sub$celltype <- tmp.sub$type
dat.res <- merge(dat,c(myeloid.sub,tcell.sub,tmp.sub))

 saveRDS(dat.res,file="/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")

###########################################################################################################################################################
# run cellchat 
# use type is three group and big celltype 
# celltype is small celltype and type_group is cell type 
#==========================================================================================================================================================
library(Seurat)
library(CellChat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
data.input <- GetAssayData(dat, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(dat@meta.data) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
levels(cellchat@idents)

CellChatDB <- CellChatDB.human
# use all DB 
CellChatDB.use <- CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel

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

saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subTumor_celltype.cellchat.RDS")











# 2021-7-12 
# check result 
##############################################################################################################################################
library(Seurat)
library(CellChat)

cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subTumor_type.cellchat.RDS")
groupSize <- table(cellchat@idents)
netVisual_circle(cellchat@net$count, sources.use=4,targets.use = c(1,2,6,7,8),vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

tmp <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1,2,7), remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
#rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))







#======================================================================================================================================
# check Treg communication
# 2021-7-15
#======================================================================================================================================
library(Seurat)
library(CellChat)

cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subTumor_celltype.cellchat.RDS")
groupSize <- table(cellchat@idents)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

tmp <- netVisual_bubble(cellchat, sources.use =c(7,9), targets.use = 8, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
#rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))














######################################################################################################################################################################
# 2021-7-29 
# use gene expression specifity to check interaction 
#=====================================================================================================================================================================
library(Seurat)
source("/public/workspace/lily/software/specInt.R")
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
subdat <- subset(dat,cells=which(dat$celltype%in%c("group09","group16","group47","Oligodendrocyte","Endothelial","Fibroblast")))
# use Function 
rs = specInt(as.matrix(subdat[["RNA"]]@data),subdat$celltype,celldb[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')],
    c("group16"),c("Fibroblast")
)
rs$pairname <- celldb[rownames(rs),"Pair.Name"]
rs.f = rs[which(rs$L.mean>0.1 & rs$R.mean>0.1),]
rs.f$L.fc = rs.f$L.mean/rs.f$nL.mean
rs.f$R.fc = rs.f$R.mean/rs.f$nR.mean
head(rs.f[order(apply(rs.f[,c('L.fc','R.fc')],1,min),decreasing=T),],10)











#========================================================================================================================================
# 2021-8-18
# 2021-9-8
library(Seurat)
# source("/public/workspace/lily/software/specInt.R")
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
# subdat <- subset(dat,cells=which(dat$celltype%in%c("group09","group16","group47","MDM","Treg")))
# # use Function 
# rs = specInt(as.matrix(subdat[["RNA"]]@data),subdat$celltype,celldb[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')],
#     c("group09"),c("MDM")
# )

source("/public/workspace/lily/software/specInt.R")
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
subdat <- subset(dat,cells=which(dat$celltype%in%c("group09","group16","group47","Treg")))
genepair <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
genepair$Pair.Name <- paste0(genepair$ligand,"_",genepair$receptor)
rs = specInt(as.matrix(subdat[["RNA"]]@data),subdat$celltype,genepair,
    c("group09"),c("Treg")
)

rs$pairname <- genepair[rownames(rs),"Pair.Name"]
rs.f = rs[which(rs$L.mean>0.1 & rs$R.mean>0.1),]
rs.f$L.fc = rs.f$L.mean/rs.f$nL.mean
rs.f$R.fc = rs.f$R.mean/rs.f$nR.mean
head(rs.f[order(apply(rs.f[,c('L.fc','R.fc')],1,min),decreasing=T),],10)
















############################################################################################################################################
# 2021-9-13
# analysis macrophage factor 
#===========================================================================================================================================
library(Seurat)

tumor.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
tumor.e@active.ident <- factor(tumor.e$tmp.group)

tumor.d <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
tumor.d@active.ident <- factor(tumor.d$tmp.group)

# LX701 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor_LX701.RDS")
# DefaultAssay(LX701) <- "RNA"
# LX701@active.ident <- factor(LX701$tmp.group)

# LX255B <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor_LX255B.RDS")
# DefaultAssay(LX255B) <- "RNA"
# LX255B@active.ident <- factor(LX255B$tmp.group)

# LX681 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor_LX681.RDS")
# DefaultAssay(LX681) <- "RNA"
# LX681@active.ident <- factor(LX681$tmp.group)

lcbm <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(lcbm) <- "RNA"
lcbm@active.ident <- factor(lcbm$group)

gse123902 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
DefaultAssay(gse123902) <- "RNA"
gse123902@active.ident <- factor(gse123902$tmp.group)

features <- c("CSF1","CCL2","CCL20","CCL5","CCL3","ICAM1")
# do not have CXCL12 and LGALS3 not each sample have this gene expression 

features <- c("CD274","LGALS9","PVR","VTCN1","IDO1")
tmp.e <- AverageExpression(tumor.e,assay="RNA",features=features)$RNA
tmp.d <- AverageExpression(tumor.d,assay="RNA",features=features)$RNA
tmp.lx701 <- AverageExpression(LX701,assay="RNA",features=features)$RNA
tmp.lx255b <- AverageExpression(LX255B,assay="RNA",features=features)$RNA
tmp.lx681 <- AverageExpression(LX681,assay="RNA",features=features)$RNA
tmp.lcbm <- AverageExpression(lcbm,assay="RNA",features=features)$RNA



# remove unclassify and do scale
tmp.e <- tmp.e[,-4]
tmp.d <- tmp.d[,-4]
tmp.lx701 <- tmp.lx701[,-4]
tmp.lx255b <- tmp.lx255b[,-4]
tmp.lx681 <- tmp.lx681[,-4]

exp.e <- t(apply(tmp.e,1,function(x){scale(as.numeric(x))}))
colnames(exp.e) <- colnames(tmp.e)
exp.d <- t(apply(tmp.d,1,function(x){scale(as.numeric(x))}))
colnames(exp.d) <- colnames(tmp.d)
exp.lx701 <- t(apply(tmp.lx701,1,function(x){scale(as.numeric(x))}))
colnames(exp.lx701) <- colnames(tmp.lx701)
exp.lx255b <- t(apply(tmp.lx255b,1,function(x){scale(as.numeric(x))}))
colnames(exp.lx255b) <- colnames(tmp.lx255b)
exp.lx681 <- t(apply(tmp.lx681,1,function(x){scale(as.numeric(x))}))
colnames(exp.lx681) <- colnames(tmp.lx681)
exp.lcbm <- t(apply(tmp.lcbm,1,function(x){scale(as.numeric(x))}))
colnames(exp.lcbm) <- colnames(tmp.lcbm)


#====================================================================================================================================
# plot heatmap for each sample 
library(pheatmap)

# pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_recruit_mac.pdf",useDingbats=F)
# pheatmap(tmp.e[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.e")
# pheatmap(tmp.d[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.d")
# pheatmap(tmp.lx701[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.lx701")
# pheatmap(tmp.lx255b[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.255b")
# pheatmap(tmp.lx681[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.lx681")
# pheatmap(tmp.lcbm,scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,main="tmp.lcbm")
# dev.off()

#gse123902 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
# add percentage info 
percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
            tmpp <- c(0,0,0,0)
            names(tmpp) <- names(table(dat@meta.data[,group]))
           res.list[[i]] <- tmpp
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}

res.e.list <- percent_feature(tumor.e,features,group="tmp.group")
res.d.list <- percent_feature(tumor.d,features,group="tmp.group")
# res.lx681.list <- percent_feature(LX681,features,group="tmp.group")
# res.lx701.list <- percent_feature(LX701,features,group="tmp.group")
# res.lx255b.list <- percent_feature(LX255B,features,group="tmp.group")
res.lcbm.list <- percent_feature(lcbm,features,group="group")
res.gse123902.list <- percent_feature(gse123902,features,group="tmp.group")






#===================================================================================================================================
# just use percentage to plot result 
# 2021-9-13
#===================================================================================================================================
library(Seurat)
list2mat <- function(res.list){
    res.dat <- matrix(unlist(res.list),ncol=length(res.list[[1]]),byrow=T)
    rownames(res.dat) <- unique(sapply(strsplit(as.vector(names(unlist(res.list))),"\\."),function(x){x[[1]]}))
    colnames(res.dat) <- unique(sapply(strsplit(as.vector(names(unlist(res.list))),"\\."),function(x){x[[2]]}))
    return(res.dat)
}
mat.e <- list2mat(res.e.list)
mat.d <- list2mat(res.d.list)
# mat.lx701 <- list2mat(res.lx701.list)
# mat.lx681 <- list2mat(res.lx681.list)
# mat.lx255b <- list2mat(res.lx255b.list)
mat.lcbm <- list2mat(res.lcbm.list)
mat.gse123902 <- list2mat(res.gse123902.list)

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_recruit_mac_percent_heatmap.pdf",useDingbats=F)
pheatmap(mat.e[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.e")
pheatmap(mat.d[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.d")
# pheatmap(mat.lx701[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
# pheatmap(mat.lx681[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
# pheatmap(mat.lx255b[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
pheatmap(mat.lcbm[,c(3,1,2)],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.lcbm")
pheatmap(mat.gse123902[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.gse123902")
dev.off()



















# combine data
combine_data <- function(res.list,tmp.dat){
    library(reshape2) 
    tmp.dat <- data.frame(tmp.dat)
    res.dat <- data.frame(percent=unname(unlist(res.list)),name=names(unlist(res.list)))
    tmp.dat$gene <- rownames(tmp.dat)
    tmp.res.dat <- melt(tmp.dat,id="gene")
    tmp.res.dat$name <- paste0(tmp.res.dat$gene,".",tmp.res.dat$variable)
    colnames(tmp.res.dat)[3] <- ("Exp")
    res.e <- merge(res.dat,tmp.res.dat[,c("Exp","name")])

    # split
    res.e$gene <- sapply(strsplit(as.vector(res.e$name),"\\."),function(x){x[[1]]})
    res.e$group <- sapply(strsplit(as.vector(res.e$name),"\\."),function(x){x[[2]]})
    return(res.e)
}

res.e <- combine_data(res.e.list,exp.e)
res.d <- combine_data(res.d.list,exp.d)
res.lx255b <- combine_data(res.lx255b.list,exp.lx255b)
res.lx701 <- combine_data(res.lx701.list,exp.lx701)
res.lx681 <- combine_data(res.lx681.list,exp.lx681)
res.lcbm <- combine_data(res.lcbm.list,exp.lcbm)








# plot result

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_recruit_mac_bubble.pdf",useDingbats=F)
ggplot(res.e,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.e")

ggplot(res.d,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") + theme_bw() + ggtitle("tmp.d") 

ggplot(res.lx701,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lx701") 

ggplot(res.lx681,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lx681")

ggplot(res.lx255b,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lx255b")

res.lcbm$group <- factor(res.lcbm$group,levels=c("group47","group09","group16"))
ggplot(res.lcbm,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lcbm")
dev.off()



library(ggpubr)
p1 <- ggplot(res.d,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") + theme_bw() + ggtitle("tmp.d") 

res.lcbm$group <- factor(res.lcbm$group,levels=c("group47","group09","group16"))
p2 <- ggplot(res.lcbm,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lcbm")

p3 <- ggplot(res.lx681,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lx681")

p4 <- ggplot(res.lx255b,aes(x=group,y=gene,size=percent,color=Exp))+geom_point() + scale_size_continuous(range = c(0,8)) +
    scale_colour_gradient2(low="steelblue",mid="white",high="#E41A1C") +theme_bw() + ggtitle("tmp.lx255b")

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_recruit_mac_bubble_select.pdf",useDingbats=F)
ggarrange(p1,p2,p3,p4,ncol = 2,nrow=2,
          common.legend = T)
dev.off()









#=======================================================================================================================================================
# 2021-9-14
# verify lung brain MDM MG signature 
#=======================================================================================================================================================
# 1. GSE14108
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

features <- c("CSF1","CCL2","CCL20","CCL5","CCL3","ICAM1")
res.list <- list()
for(i in 1:length(features)){
    if(features[i]%in%rownames(dat.BM)){
        gene_exp <- as.numeric(dat.BM[features[i],])
        tmp <- apply(mod[,6:10],2,function(x){
            c(cor.test(x,gene_exp,method="spearman")$estimate, cor.test(x,gene_exp,method="spearman")$p.value)
        })
    }else{
        tmp <- c(0,0)
    }
    res.list[[i]] <- tmp
    names(res.list)[i] <- features[i]  
}
# do not detec CCL3 ,remove
res.list$CCL3 <- NULL
tmp <- sapply(res.list,function(x){as.matrix(x)})
rownames(tmp) <- c("Lung.cor","Lung.p","Brain.cor","Brain.p","BMS.cor","BMS.p","BMDM.cor","BMDM.p","MG.cor","MG.p")
tmp <- data.frame(t(tmp))
tmp$gene <- rownames(tmp)


# tmp$color <- "NS."
# tmp$color[which(rownames(tmp)%in%c("CCL20","ICAM1"))] <- "Lung.BMS"
# tmp$color[which(rownames(tmp)%in%c("CCL3"))] <- "Brain"
# tmp$color[which(rownames(tmp)%in%c("CCL2"))] <- "Brain.BMS"
# tmp$color[which(rownames(tmp)%in%c("CCL5","CSF1"))] <- "BMS"


tmp$shape <- "NS."
tmp$shape[which(rownames(tmp)%in%c("CSF1","CCL2","CCL20","CCL5"))] <- "MDM"

library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GSE14108_factor_BMS.pdf",useDingbats=F)
ggplot(tmp,aes(x=Lung.cor,y=BMS.cor,shape=shape,size=I(5)))+ geom_point() + geom_text(aes(label=gene, y=BMS.cor+0.05)) +
    scale_shape_manual(values=c(19,15))
dev.off()




# 2. E-MTAB 
# use this because have CCL3 

dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

features <- c("CSF1","CCL2","CCL20","CCL5","CCL3","ICAM1")
res.list <- list()
for(i in 1:length(features)){
    if(features[i]%in%rownames(dat.BM)){
        gene_exp <- as.numeric(dat.BM[features[i],])
        tmp <- apply(mod[,6:10],2,function(x){
            c(cor.test(x,gene_exp,method="spearman")$estimate, cor.test(x,gene_exp,method="spearman")$p.value)
        })
    }else{
        tmp <- c(0,0)
    }
    res.list[[i]] <- tmp
    names(res.list)[i] <- features[i]  
}

# plot dot plot to show result 
tmp <- sapply(res.list,function(x){as.matrix(x)})
rownames(tmp) <- c("Lung.cor","Lung.p","Brain.cor","Brain.p","BMS.cor","BMS.p","BMDM.cor","BMDM.p","MG.cor","MG.p")
tmp <- data.frame(t(tmp))
tmp$gene <- rownames(tmp)

tmp$color <- "NS."
tmp$color[which(rownames(tmp)%in%c("CCL20","ICAM1"))] <- "Lung.BMS"
tmp$color[which(rownames(tmp)%in%c("CCL3"))] <- "Brain"
tmp$color[which(rownames(tmp)%in%c("CCL2"))] <- "Brain.BMS"
tmp$color[which(rownames(tmp)%in%c("CCL5","CSF1"))] <- "BMS"


tmp$shape <- "NS."
tmp$shape[which(rownames(tmp)%in%c("CCL2","CCL5","CCL3"))] <- "MDM.MG"
tmp$shape[which(rownames(tmp)%in%c("CSF1","ICAM1","CCL20"))] <- "MDM"


#########################################################################################################################################################
# ggplot  show result 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/E-MATB_factor_BMS.pdf",useDingbats=F)
ggplot(tmp,aes(x=Lung.cor,y=BMS.cor,shape=shape,color=color,size=I(5)))+ geom_point() + geom_text(aes(label=gene, y=BMS.cor+0.05)) +
    scale_colour_manual(values=c("#FFC719","#00A5DC","#00b2a9","#a0ac48"))
dev.off()



# try to plot heatmap
#===================================================================================================================================================== 
library(pheatmap)
dat <- dat.BM[features,]
ann <- data.frame(row.names=colnames(dat),Lungsig= mod[,6],Brainsig=mod[,7],BMS=mod[,8],BMDM=mod[,9],MG=mod[,10])
ann <- ann[order(ann$BMS,decreasing=T),]

pheatmap(dat[,rownames(ann)],annotation_col=ann,scale="column")





















































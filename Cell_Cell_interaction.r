

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
dat@active.ident <- factor(dat$celltype)
tmp <- AverageExpression(dat,assay="RNA")$RNA
tmp$M.T.FC <- tmp$MDM/tmp$tumor
tmp$T.M.FC <- tmp$tumor/tmp$MDM
# tmp <- tmp[order(tmp$M.T.FC,decreasing=T),]
tmp.f <- tmp[-which(is.infinite(tmp$T.M.FC)|is.infinite(tmp$M.T.FC)),]
MDM.h <- rownames(tmp.f)[which(tmp.f$M.T.FC>2)]
tumor.f <- rownames(tmp.f)[which(tmp.f$T.M.FC>100)]
# 2021-6-10 
library(CellChat)
CellChatDB <- CellChatDB.human

info <-  data.frame(ligand=sapply(strsplit(as.vector(CellChatDB$interaction$interaction_name_2)," - "),function(x){x[1]}),
receport=sapply(strsplit(as.vector(CellChatDB$interaction$interaction_name_2)," - "),function(x){x[2]}),stringsAsFactors=F)

for(i in 1:nrow(info)){
    if(length(grep("\\)|\\(",info[i,2]))>0){
        tmp <- gsub("\\(|\\)|^ ","",info[i,2])
        pair <- strsplit(tmp,"\\+")[[1]]
        for(j in 1:length(pair)){
            value <- c(info[i,1],pair[j])
            info <- rbind(info,value)
        }
    }
}
info.f <- info[-grep("\\+",info$receport),] # get gene pair 
colnames(info.f) <- c("ligand","receptor")
saveRDS(info.f,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
# get gene expression 

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
dat@active.ident <- factor(dat$celltype)
tmp <- AverageExpression(dat,assay="RNA")$RNA
tmp$M.T <- tmp$MDM-tmp$tumor
tmp$T.M <- tmp$tumor-tmp$MDM
tmp$gene <- rownames(tmp)

# LUAD 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")
dat@active.ident <- factor(dat$celltype)
tmp <- AverageExpression(dat,assay="RNA")$RNA
tmp$M.T <- tmp$MDM-tmp$tumor
tmp$T.M <- tmp$tumor-tmp$MDM
tmp$gene <- rownames(tmp)


# GBM 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")
dat@active.ident <- factor(dat$celltype)
tmp <- AverageExpression(dat,assay="RNA")$RNA
tmp$M.T <- tmp$MDM-tmp$GBM_tumor
tmp$T.M <- tmp$GBM_tumor-tmp$MDM
tmp$gene <- rownames(tmp)


genepair <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
# genepair$ligand <- gsub(" $","",genepair$ligand)
# genepair$receptor <- gsub("^ ","",genepair$receptor)
# saveRDS(genepair,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
tmp1 <- merge(tmp[,c("M.T","gene")],genepair,by.x="gene",by.y="ligand",all.y=T) 
colnames(tmp1)[1] <- "ligand"
tmp2 <- merge(tmp[,c("T.M","gene")],tmp1,by.x="gene",by.y="receptor",all.y=T) 
colnames(tmp2)[1] <- "receptor"
tmp2[which(tmp2$M.T>0.5&tmp2$T.M>0.5),]










###################################################################################################################
# 2021-6-11
# 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
dat@active.ident <- factor(dat$celltype)
tmp <- AverageExpression(dat,assay="RNA")$RNA
tmp$M.T <- tmp$MDM-tmp$tumor
tmp$T.M <- tmp$tumor-tmp$MDM
tmp$gene <- rownames(tmp)

dat$class1 <- as.vector(dat$celltype)
dat$class1[which(dat$celltype%in%c("MG","Monocyte","tumor"))] <- "other"
dat@active.ident <- factor(dat$class1)
tmp1 <- AverageExpression(dat,assay="RNA")$RNA
tmp1$MDM.h <- tmp1$MDM -tmp1$other

dat$class2 <- as.vector(dat$celltype)
dat$class2[which(dat$celltype%in%c("MG","Monocyte","MDM"))] <- "other"
dat@active.ident <- factor(dat$class2)
tmp2 <- AverageExpression(dat,assay="RNA")$RNA
tmp2$tumor.h <- tmp2$tumor -tmp2$other

tmp <- cbind(tmp1$MDM.h,tmp2$tumor.h)
colnames(tmp) <- c("M.T","T.M")
rownames(tmp) <- rownames(tmp1)
tmp <- data.frame(tmp)
tmp$gene <- rownames(tmp)


genepair <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
# genepair$ligand <- gsub(" $","",genepair$ligand)
# genepair$receptor <- gsub("^ ","",genepair$receptor)
# saveRDS(genepair,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/gene_pair.RDS")
tmp1 <- merge(tmp[,c("M.T","gene")],genepair,by.x="gene",by.y="receptor",all.y=T) 
colnames(tmp1)[1] <- "receptor"
tmp2 <- merge(tmp[,c("T.M","gene")],tmp1,by.x="gene",by.y="ligand",all.y=T) 
colnames(tmp2)[1] <- "ligand"

tmp2[which(tmp2$M.T>0.5&tmp2$T.M>0.5),]



tmp1 <- merge(tmp[,c("M.T","gene")],genepair,by.x="gene",by.y="ligand",all.y=T) 
colnames(tmp1)[1] <- "ligand"
tmp2 <- merge(tmp[,c("T.M","gene")],tmp1,by.x="gene",by.y="receptor",all.y=T) 
colnames(tmp2)[1] <- "receptor"
tmp2[which(tmp2$M.T>0.5&tmp2$T.M>0.5),]


# get gene pair names

tmp.log <- log2(tmp[,c("MDM","tumor")]+1)
tmp.log$group <- "no"
tmp.log$group[which(rownames(tmp.log)%in%c("SDC4","OSMR","EGFR","ITGA3","IL1R1","SDC1"))] <- "receptor"
tmp.log$group[which(rownames(tmp.log)%in%c("FN1","THBS1","OSM","IL1B","HBEGF"))] <- "ligand"
tmp.log$group <- factor(tmp.log$group,levels=c("receptor","ligand","no"))

tmp.log$gene <- ""
tmp.log$gene[which(rownames(tmp.log)%in%c("SDC4","OSMR","EGFR","ITGA3","IL1R1","SDC1"))] <- rownames(tmp.log)[which(rownames(tmp.log)%in%c("SDC4","OSMR","EGFR","ITGA3","IL1R1","SDC1"))]
tmp.log$gene[which(rownames(tmp.log)%in%c("FN1","THBS1","OSM","IL1B","HBEGF"))] <- rownames(tmp.log)[which(rownames(tmp.log)%in%c("FN1","THBS1","OSM","IL1B","HBEGF"))]


pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/MDM2Tumor.interaction.pdf",useDingbats=F)
ggplot(tmp.log,aes(x=MDM,y=tumor,colour=group,group=group)) + geom_point() + theme_classic() +
    scale_colour_manual(values=c("#007cc0","#ffb310","grey")) +
    geom_text(aes(label = gene), size = 3)    

dev.off()










#######################################################################################################################################
# 2021-6-17
# CellChat result check 
#======================================================================================================================================
# 1. MDM2Tumor 
library(Seurat)
library(CellChat)


cellchat1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.CellChat.RDS")
cellchat2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/GBM.CellChat.RDS")
cellchat3 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LUAD.CellChat.RDS")

# LCBM 
tmp <- netVisual_bubble(cellchat1, sources.use = 3, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "LCBM"
tmp1 <- tmp.res
# LUAD
tmp <- netVisual_bubble(cellchat3, sources.use = 1, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "LUAD"
tmp3 <- tmp.res
# GBM 
tmp <- netVisual_bubble(cellchat2, sources.use = 2, targets.use = 1, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "GBM"
tmp2 <- tmp.res

# merge interaction data 
tmp.f <- merge(tmp1,tmp2,all.x=T,all.y=T,by="interaction_name_2")
colnames(tmp.f)[2:5] <- c("prob.LCBM","group.LCBM","prob.GBM","group.GBM")
tmp.rs <- merge(tmp.f,tmp3,all.x=T,all.y=T,by="interaction_name_2")

sub.rs <- tmp.rs[,c("interaction_name_2","prob.LCBM","prob.GBM","prob")]
colnames(sub.rs)[4] <- "prob.LUAD"



tmp.lung <- sub.rs[which(sub.rs$prob.LUAD>0&sub.rs$prob.LCBM>0&is.na(sub.rs$prob.GBM)),]
rownames(tmp.lung) <- tmp.lung$interaction_name_2
tmp.lung <- tmp.lung[order(tmp.lung$prob.LCBM,decreasing=T),]
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM.LUAD.CellChat.MDM.tumor.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.lung[,2:4],na_col="grey",cluster_rows=F,cluster_cols=F)
dev.off()

tmp.brain <- sub.rs[which(sub.rs$prob.GBM>0&sub.rs$prob.LCBM>0&is.na(sub.rs$prob.LUAD)),]
rownames(tmp.brain) <- tmp.brain$interaction_name_2
tmp.brain <- tmp.brain[order(tmp.brain$prob.LCBM,decreasing=T),]
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM.GBM.CellChat.MDM.tumor.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.brain[,2:4],na_col="grey",cluster_rows=F,cluster_cols=F)
dev.off()







# 2021-7-6
# 2. Tumor2MDM 
#####################################################################################################################################
library(Seurat)
library(CellChat)


cellchat1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.CellChat.RDS")
cellchat2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/GBM.CellChat.RDS")
cellchat3 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LUAD.CellChat.RDS")

# LCBM 
tmp <- netVisual_bubble(cellchat1, sources.use = 4, targets.use = 1, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "LCBM"
tmp1 <- tmp.res
# LUAD
tmp <- netVisual_bubble(cellchat3, sources.use = 4, targets.use = 1, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "LUAD"
tmp3 <- tmp.res
# GBM 
tmp <- netVisual_bubble(cellchat2, sources.use = 1, targets.use = 2, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "GBM"
tmp2 <- tmp.res

# merge interaction data 
tmp.f <- merge(tmp1,tmp2,all.x=T,all.y=T,by="interaction_name_2")
colnames(tmp.f)[2:5] <- c("prob.LCBM","group.LCBM","prob.GBM","group.GBM")
tmp.rs <- merge(tmp.f,tmp3,all.x=T,all.y=T,by="interaction_name_2")

sub.rs <- tmp.rs[,c("interaction_name_2","prob.LCBM","prob.GBM","prob")]
colnames(sub.rs)[4] <- "prob.LUAD"


tmp.lung <- sub.rs[which(sub.rs$prob.LUAD>0&sub.rs$prob.LCBM>0&is.na(sub.rs$prob.GBM)),]
rownames(tmp.lung) <- tmp.lung$interaction_name_2
tmp.lung <- tmp.lung[order(tmp.lung$prob.LCBM,decreasing=T),]
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM.LUAD.CellChat.Tumor.mdm.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.lung[,2:4],na_col="grey",cluster_rows=F,cluster_cols=F)
dev.off()

tmp.brain <- sub.rs[which(sub.rs$prob.GBM>0&sub.rs$prob.LCBM>0&is.na(sub.rs$prob.LUAD)),]
rownames(tmp.brain) <- tmp.brain$interaction_name_2
tmp.brain <- tmp.brain[order(tmp.brain$prob.LCBM,decreasing=T),]
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM.GBM.CellChat.Tumor.mdm.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.brain[,2:4],na_col="grey",cluster_rows=F,cluster_cols=F)
dev.off()













# 2021-7-6
# calculate LCBM samples Tumor with Myeloid cell cell communication 
####################################################################################################################################
library(Seurat)
library(CellChat)
# make a data 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))
tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tumor$celltype <- tumor$group
tmp.res <- merge(sub.dat,tumor)
sub.dat.f <- subset(tmp.res,cells=which(tmp.res$celltype%in%c("MDM","MG","Monocyte","group09","group16","group47")))
saveRDS(sub.dat.f,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.RDS")

# by the way ,do cell chat analysis 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.RDS")
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

saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.cellchat.RDS")





#===================================================================================================================================
# check result 
# 0. interaction number differneces 
library(CellChat)
cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.cellchat.RDS")
groupSize <- table(cellchat@idents)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_bubble(cellchat, sources.use = 1, targets.use = 2, remove.isolate = FALSE,return.data=F)


cellchat@net$count -> tmp
tmp[1:3,1:3]
tmp[1:3,1:3] <- 0
tmp[4:6,4:6] <- 0

rownames(tmp) <- paste0(rownames(tmp),".source")
colnames(tmp) <- paste0(colnames(tmp),".target")
groups <- rep(c("group09","group16","group47","MDM","MG","Monocyte"),2)
names(groups) <- c(rownames(tmp),colnames(tmp))
cols <- c("#fbbc05","#fbbc05","#34a853","#34a853","#ea4335","#ea4335")
names(cols) <-c("Myeloid.source","Myeloid.target","T_cell.source","T_cell.target","Tumor_cell.source","Tumor_cell.target")

netVisual_chord_cell_internal(tmp,color.use=cols,sources.use=rownames(tmp),targets.use=colnames(tmp),group=groups,big.gap=20,small.gap=0)



# 1. tumor interaction with Myeloids
library(CellChat)
cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.cellchat.RDS")
groupSize <- table(cellchat@idents)
tmp <- netVisual_bubble(cellchat, sources.use = 4, targets.use = 1:3, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))





















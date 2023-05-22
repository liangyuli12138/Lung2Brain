
# 2021-6-23 
# this program is used to analysis T cell 
# result plot should be in /public/workspace/lily/Lung2Brain/Version5/T_cell/plot
# 1. T cell TSNE 
# 2. bar plot show cell percentage change in different group
# 3. T cell exhausted signature calculate ?





#1. T cell subgroup 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
#        Bcell               CD4 naive           CD8 cytotoxic
#                    1189                    2744                    2697
# CD8 exhausted/cytotoxic                  NKcell                  Plasma
#                    3035                     702                     680
#                     Tfh                    Treg                Undefine
#                     581                    1567                     609

cols <- c("#54B0E4","#4DAF4A","#F29403","#c6af92","#B3DE69","#00CDD1","#BC9DCC","#5E4FA2","#999999")
DimPlot(dat,group.by="celltype",cols=cols)
pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/Tcell_tsne.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype",cols=cols)
dev.off()



# 2.cell percentage 
library(reshape)
library(ggplot2)
library(ggalluvial)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
tmp <- table(dat$type_group,dat$celltype)
tmp.f <- tmp[,-grep("Undefine",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,level=c("GBM","LCBM","LC"))
pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/Tcell_celltype.pdf",useDingbats=F)
cols <- c("#54B0E4","#4DAF4A","#F29403","#c6af92","#B3DE69","#00CDD1","#BC9DCC","#5E4FA2")
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()







# 3. TCGA and GSE161116 to verify 
# GSE161116 show converse trend (brain metastasis have less Treg expression [FOXP3,IL2RA])
# so maybe use GSE126548 to verify 
#========================================================================================================
# Treg signature 
# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2006-7-7-r54/figures/2
# gene[-which(gene%in%rownames(dat))] # change some gene alias 
# [1] "BHLHB2"   "CEB1"     "GPR2"     "HLA-DRB3" "G1P2"
# BHLHE40 HERC5 CCR10 no ISG15
# gene <- c("FOXP3","SDC4","NINJ2",'PTTG1','TIAF1','TRIB1','S100A10','GBP2','GATA3','IL2RA',
# 'BHLHE40','HERC5','CTLA4','TFRC','HLA-DMA','AKAP2','TNFRSF1B','CCR5','CCR10','IL2RB',
# 'SHMT2','HLA-DRB1','HLA-DRB3','TP53INP1','GBP5','EPSTI1','LGALS3','SLAMF1','TRAF1',
# 'LGALS1','S100A4',"ISG15")

# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(gene,"Treg",out="/public/workspace/lily/MOD_file/Treg.mod") # make a mod file 



# TCGA data calculate BMS and Treg signature 
dat <-readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

# plot result 
pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/TCGA_BMS_Treg.pdf",useDingbats=F)
plot(mod$BMS_update_norm,mod$Treg_norm)
abline(lm(mod$BMS_update_norm~mod$Treg_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.37))
dev.off()





# GSE161116
# GSE161116 show not significant 
#===========================================================================================================================================
# dat <- read.table("~/metastasis/data/verify/GSE161116/GSE161116_series_matrix.txt",sep="\t",header=T,comment.char="!")
# rownames(dat) <- dat$ID_REF
# dat$ID_REF <- NULL
# info <- as.vector(read.table("~/metastasis/data/verify/GSE161116/sampleinfo.txt"))[-1]
# names(info) <- NULL

# ann <- data.frame(sample=colnames(dat),patient=sapply(sapply(info,function(x){strsplit(as.vector(x)," ")}),function(y){paste0(y[1],y[2])}),
#     group=sapply(sapply(info,function(x){strsplit(as.vector(x)," ")}),function(y){paste0(y[3])})
# )

# # calculate BMS 
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)


# GSE126548 to verify result
# show correlation positive but not significant 
#===========================================================================================================================================
# load("~/metastasis/data/verify/GSE126548/GSE126548_rpkm.RData")
# dat <- rpkm
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)














# 4. Tumor and T cell analysis 
# use CellChat to analysis Tumor cells and Tcells 
#===============================================================================================================================
library(Seurat)
library(CellChat)

tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tumor$celltype <- tumor$group
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
tcell <- subset(tmp,cells=which(tmp$celltype%in%c("CD4 naive","CD8 cytotoxic","CD8 exhausted/cytotoxic","NKcell","Tfh","Treg")&tmp$type_group=="LCBM"))
dat <- merge(tumor,tcell)

# set a CellChat object 
# Run CellChat 
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

saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/T_cell/LCBM_tumor_tcell.cellchat.RDS")

#===================================================================================================









################################################################################################################################
# 5. another way to analysis T cell and Tumor cell
# 2021-6-25
# dat <- read.table("/public/workspace/lily/Lung2Brain/Version5/PairsLigRec.txt",sep="\t",header=T)
library(Seurat)

tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tumor$celltype <- tumor$group
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
tcell <- subset(tmp,cells=which(tmp$celltype%in%c("CD4 naive","CD8 cytotoxic","CD8 exhausted/cytotoxic","NKcell","Tfh","Treg")&tmp$type_group=="LCBM"))
dat <- merge(tumor,tcell)

gene <- list()
group <- names(table(dat$celltype))
for(i in 1: length(group)){
    tmp <- NULL
    if(i %in% c(1:3,7:9)){
        tmp <- FindMarkers(tcell,assay="RNA",ident.1=group[i],min.pct=0.1,only.pos=T,group.by="celltype",logfc.threshold=0.1)
    }

    if(i %in% c(4:6)){
        tmp <- FindMarkers(tumor,assay="RNA",ident.1=group[i],min.pct=0.1,only.pos=T,group.by="celltype",logfc.threshold=0.1)
    }
    
    gene[[i]] <- tmp
    names(gene)[i] <- gsub(" |/",".",group[i])

}

ann <- read.table("/public/workspace/lily/Lung2Brain/Version5/PairsLigRec.txt",sep="\t",header=T)


ann.f <- ann[which(ann$Receptor.ApprovedSymbol%in%rownames(gene$Treg)&ann$Ligand.ApprovedSymbol%in%rownames(gene$group09)),]
ann.f[,c(1,2,4)]



ann.f <- ann[which(ann$Receptor.ApprovedSymbol%in%rownames(gene$Treg)&ann$Ligand.ApprovedSymbol%in%rownames(gene$CD8.cytotoxic)),]












###############################################################################################################################################
# plot T cell results
# 1. scatter plot show correlation in GSE14108 Treg and BMS update 
# 2. scatter plot show correlation in GSE14108 Trrg and Lung_gene signature 
# 3. circle plot show subgroup interation number with Treg 
# 4. heatmap show in group47 and grou09 interaction with
# 5. correlation for LGAS9 and Treg
#=================================================================================================================================
# 2021-7-7
# 1 and 2 GSE 14108 data to verify 
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Treg","Lung_gene","Brain_gene","BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE14108_LCBM_treg_lung_BMS.pdf",useDingbats=F)
plot(mod$BMS_update_norm,mod$Treg_norm,mian="BMS with Treg")
abline(lm(mod$BMS_update_norm~mod$Treg_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.455," pvalue=",0.016))
#############################################################
plot(mod$Lung_gene_norm,mod$Treg_norm,mian="Lung with Treg")
abline(lm(mod$Lung_gene_norm~mod$Treg_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.574," pvalue=",0.001))

dev.off()




# 3. and 4. CellChat result show circle plot 
library(Seurat)
library(Cellchat)

# analysis result and check 
# 2021-7-7
library(Seurat)
library(CellChat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/LCBM_tumor_tcell.cellchat.RDS")
groupSize <- table(dat@idents)

# 3. circle plot show cell chat interaction 
netVisual_circle(dat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

tmp <- dat@net$count
tmp.f <- tmp[c(4,6,9),c(4,6,9)]
tmp.f[c(1:2),c(1:2)] <- 0
tmp.f[3,3] <- 0

rownames(tmp.f) <- paste0(rownames(tmp.f),".source")
colnames(tmp.f) <- paste0(colnames(tmp.f),".target")
groups <- rep(c("group09","group47","Treg"),2)
names(groups) <- c(rownames(tmp.f),colnames(tmp.f))


netVisual_chord_cell_internal(tmp.f,sources.use=rownames(tmp.f),targets.use=colnames(tmp.f),group=groups,big.gap=20,small.gap=0)




# 4. pheatmap show result 
#===================================================================================================================================
# netVisual_bubble(dat, sources.use = c(4,5,6), targets.use = c(9), remove.isolate = FALSE)
tmp <- netVisual_bubble(dat, sources.use = c(4,5,6), targets.use = c(9), remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
rownames(rs) <- rs$source.target
rs$source.target <- NULL
rs <- t(rs)
rs <- data.frame(rs)

tmp09 <- rs[which(rs$group09....Treg>0&is.na(rs$group47....Treg)),]
tmp09 <- tmp09[order(tmp09$group09....Treg,decreasing=T),]

tmp47 <- rs[which(rs$group47....Treg>0&is.na(rs$group09....Treg)),]
tmp47 <- tmp47[order(tmp47$group47....Treg,decreasing=T),]


pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/LCBM_subTumor_Treg_heatmap.pdf",useDingbats=F)
pheatmap::pheatmap(tmp47,na_col="grey",cluster_rows=F,cluster_cols=F)
pheatmap::pheatmap(tmp09,na_col="grey",cluster_rows=F,cluster_cols=F)
dev.off()



# 4.1 add some data verify 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/LCBM_subTumor_LGALS9.pdf",useDingbats=F)
barplot(c(0.2554715,0.2045145,0.08139221),ylim=c(0,0.3),names=c("group09","group16","group47"),main="LGALS9 expression")
dev.off()







# 5. LGALS9 correlation with Treg signature 
# gene show no significant with Gene
# however gene correlation with gene show signifiacant
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Treg","Lung_gene","Brain_gene","BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$LGALS9 <- as.numeric(as.vector(dat.BM["LGALS9",]))

# purity 
scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GSE14108_estimate_score.gct",skip = 2,header = T)
mod$purity <- as.numeric(as.vector(scores[4,-c(1,2)]))



pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE14108_LCBM_FOXP3.pdf",useDingbats=F)
plot(mod$LGALS9,mod$FOXP3,mian="LGALS9 with FOXP3")
abline(lm(mod$FOXP3~mod$LGALS9),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.39," pvalue=",0.042))


plot(mod$LIF,mod$FOXP3,mian="LGALS9 with FOXP3")
abline(lm(mod$FOXP3~mod$LIF),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.53," pvalue=",0.004))

dev.off()






#==================================================================================================================================================
# 2021-8-15
# treg verify in TCGA and BM data (GSE14108)
# https://www.jianshu.com/p/0baac4c52ac8
#==================================================================================================================================================
source("~/software/Cibersort_R.R")
result1 <- CIBERSORT('~/software/LM22.txt','~/software/GSE14108_exp.txt', perm = 100, QN = T)
result2 <- CIBERSORT('~/software/LM22.txt','~/metastasis/data/verify/TCGA_LUAD/LUAD_RNAseq_Exp.txt', perm = 100, QN = T)
saveRDS(result2,file="~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_cibersort.RDS")

# cibersort result 
# 2021-8-18 use cibersortx result 
lcbm <- read.table("~/tmp/CIBERSORTx_GSE14108_Results.txt",sep="\t",header=T)
luad <- read.table("~/tmp/CIBERSORTx_TCGA_LUAD_Results.txt",header=T,sep="\t")
pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE14108_TCGA_treg.cibersort.pdf",useDingbats=F)
boxplot(lcbm[,10],luad[,10],outline=F,names=c("LCBM","TCGA_LUAD"))
legend("topright",legend=paste0(" pvalue<",0.001))
dev.off()



# Treg signature 
luad <- readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.luad <- mod.analyze2(as.matrix(luad),c("Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod.luad <- as.data.frame(mod.luad)

lcbm <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.lcbm <- mod.analyze2(as.matrix(lcbm),c("Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod.lcbm <- as.data.frame(mod.lcbm)


dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.mtab <- mod.analyze2(as.matrix(dat.BM),c("Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod.mtab <- as.data.frame(mod.mtab)

# cibersort result 

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE14108_TCGA_treg.sig.pdf",useDingbats=F)
boxplot(mod.lcbm[,2],mod.luad[,2],outline=F,names=c("LCBM","TCGA_LUAD"))
legend("topright",legend=paste0(" pvalue= ",0.076))
dev.off()





























#############################################################################################################################################################
# 2021-9-14
# analysis T cells cytokines
# 2021-9-23
# 2021-12-9 use BMS.group new classify  
# analysis 
#============================================================================================================================================================
library(Seurat)
lcbm <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
DefaultAssay(lcbm) <- "RNA"
lcbm.sub <- subset(lcbm,cells=which(lcbm$BMS.group%in%c("high","low")))
lcbm.sub$cell.type <- lcbm.sub$BMS.group

# myeloid
mye <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
DefaultAssay(mye) <- "RNA"
mye.lcbm <- subset(mye,cells=which(mye$type_group=="LCBM"&mye$celltype%in%c("MDM","MG","Monocyte")))
mye.lcbm$cell.type <- mye.lcbm$celltype

# lymphy
lym <- readRDS("/public/workspace/lily/Lung2Brain/Version5/T_cell/inte9_lympho.RDS")
DefaultAssay(lym) <- "RNA"
lym.lcbm <- subset(lym,cells=which(lym$type_group=="LCBM"&lym$celltype%in%c("Bcell","NKcell","Plasma")))
lym.lcbm$cell.type <- lym.lcbm$celltype

# MES
tmp <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
DefaultAssay(tmp) <- "RNA"
tmp.lcbm <- subset(tmp,cells=which(tmp$type_group=="LCBM"&tmp$type%in%c("Endothelial","Fibroblast","Oligodendrocyte")))
tmp.lcbm$cell.type <- tmp.lcbm$type

dat.res <- merge(x=lcbm.sub,y=c(mye.lcbm,lym.lcbm,tmp.lcbm))

# also analysis MDM and MG 
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/SubTumor/LCBM_subtumor_myeloid.RDS")
# tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS") # use all cell type # 2021-12-9 think this is not ok
# dat <- subset(tmp,cells=which(tmp$celltype%in%c("Endothelial","Bcell","Fibroblast","MDM","MG","Monocyte","NKcell",
#     "group09","group16","group47","Oligodendrocyte","Plasma")))
# DefaultAssay(dat) <- "RNA"
# dat@active.ident <- factor(dat$celltype)

dat <- dat.res
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$cell.type)

features <- c("CD274","LGALS9","PVR","IDO1","CCL28","CCL5","CCL8","CCL22","CCL20",
"IL6","CXCL10","CXCL9","CCL3","IFNG","CXCL11","CCL4","CCL11","IL22","CCL17","IL4","IL17F","IL5",
"IL13","IL2","IL17A","IL21","CCL2","CXCL5","CXCL1")

features <- c("LGALS9","IDO1","CCL28","CCL8","CCL22","CCL20",
"IL6","CXCL10","CXCL9","IFNG","CXCL11","IL22","IL4","IL17F","IL5",
"IL13","IL2","IL17A","IL21","CCL2","CXCL5","CXCL1")

features <- c("IL2","IFNG",
"CXCL9","CXCL10","CXCL11",
"CCL22","CCL20","CCL28",
"IDO1","LGALS9",
"CXCL1","IL6","CCL8","CCL2"
)




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

res.dat.list <- percent_feature(dat,features,group="cell.type")

list2mat <- function(res.list){
    res.dat <- matrix(unlist(res.list),ncol=length(res.list[[1]]),byrow=T)
    rownames(res.dat) <- names(res.list)
    colnames(res.dat) <- names(res.list[[1]])
    return(res.dat)
}

mat.lcbm <- list2mat(res.dat.list)
# mat.f <- mat.lcbm[-which(rowSums(mat.lcbm)==0),]
# # do some prepare work 
# mat.f[which(mat.f[]<0.01)] <- 0
# mat <- mat.f[-which(rowSums(mat.f)==0),] 

mat <- mat.lcbm[,c("NKcell","Bcell","Plasma","Fibroblast","Endothelial","Oligodendrocyte","high","low","MDM","MG","Monocyte")]

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/LCBM_Myeloid_cytokine.pdf",useDingbats=F)
library(pheatmap)
pheatmap(mat,scale="row",cluster_rows=F,cluster_cols=F,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),gaps_col=c(3,6,8))
dev.off()

# library(pheatmap)
# pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/LCBM_Myeloid_cytokine.pdf",useDingbats=F)
# pheatmap(mat.lcbm[,1:5],scale="row",cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100))
# dev.off()


# use other data to verify 
# 1. GSE14108
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
#dat.BM <- dat.BM[,20:28]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("BMS_update","BMDM_marker","MG_marker","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

# features <- c("CD274","LGALS9","PVR","VTCN1","IDO1","IDO2","CCL28","CCL5","CCL8","CCL22","CCL20")

features <- c("IL2","IFNG",
"CXCL9","CXCL10","CXCL11",
"CCL22","CCL20",
"IDO1","LGALS9",
"CXCL1","IL6","CCL8","CCL2","FOXP3","CCL28","CEBPD"
)



res.list <- list()
for(i in 1:length(features)){
    if(features[i]%in%rownames(dat.BM)){
        gene_exp <- as.numeric(dat.BM[features[i],])
        tmp <- apply(mod[,5:8],2,function(x){
            c(cor.test(x,gene_exp,method="spearman")$estimate, cor.test(x,gene_exp,method="spearman")$p.value)
        })
    }else{
        tmp <- c(0,0)
    }
    res.list[[i]] <- tmp
    names(res.list)[i] <- features[i]  
}


# 2. E-MTAB 
# use this because have CCL3 

dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("BMS_update","BMDM_marker","MG_marker","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
# features <- c("CD274","LGALS9","PVR","VTCN1","IDO1","IDO2","CCL28","CCL5","CCL8","CCL22","CCL20")
features <- c("IL2","IFNG","CEBPD",
"CXCL9","CXCL10","CXCL11",
"CCL22","CCL20",
"IDO1","LGALS9",
"CXCL1","IL6","CCL8","CCL2","FOXP3","PDCD1","CTLA4","LAG3","CCL28"
)

res.list <- list()
for(i in 1:length(features)){
    if(features[i]%in%rownames(dat.BM)){
        gene_exp <- as.numeric(dat.BM[features[i],])
        tmp <- apply(mod[,5:8],2,function(x){
            c(cor.test(x,gene_exp,method="spearman")$estimate, cor.test(x,gene_exp,method="spearman")$p.value)
        })
    }else{
        tmp <- c(0,0)
    }
    res.list[[i]] <- tmp
    names(res.list)[i] <- features[i]  
}







#===================================================================================================================================================
# 2021-9-24
# use immunity Treg signature
# seem not every good 
#===================================================================================================================================================
# tmp <- as.vector(read.table("/public/workspace/lily/tmp/tumor.Treg")[,1])
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(tmp,"tumor.Treg",out="/public/workspace/lily/MOD_file/tumor.Treg.mod") # make a mod file 



# 2021-9-24
# plot CCL20 result 
LUAD <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
luad_ann <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RDS")

rownames(luad_ann) <- gsub("-",".",rownames(luad_ann))
luad_ann <- luad_ann[colnames(LUAD),]

# LUAD three 
luad_ann$type <- rep("Unknow",nrow(luad_ann))
luad_ann[which(LUAD["CCL20",] < as.numeric(quantile(LUAD["CCL20",],0.33))),"type"] <- "CCL20.Low"
luad_ann[which(LUAD["CCL20",] > as.numeric(quantile(LUAD["CCL20",],0.67))),"type"] <- "CCL20.High"

library(ggplot2)
library(survminer)
library(survival)


pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/CCL20_surv.pdf",useDingbats=F)
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survminer::surv_fit(surv~type,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" Overall survival (months)")+labs(title="LUAD_sig_three")


surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survminer::surv_fit(surv~type,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" RF survival (months)")+labs(title="LUAD_sig_three")

dev.off()

# calculate P-value 
pairwise_survdiff(Surv(OS.time,OS)~type,data=luad_ann) # 0.043
pairwise_survdiff(Surv(RFS_time,RFS)~type,data=luad_ann) # 0.061




# CCL20 expression with Lung_signature 
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","BMDM_marker","MG_marker","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$CCL20 <- as.numeric(dat.BM["CCL20",])

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/E_MTAB_LCBM_CCL20.pdf",useDingbats=F)
plot(mod$CCL20,mod$Lung_gene_norm,mian="CCL20 with Lung sig.")
abline(lm(mod$Lung_gene_norm~mod$CCL20),col="red")
# cor.test(mod$CCL20,mod$Lung_gene_norm,method="spearman")
legend("topright",legend=paste0("rho=",0.30," pvalue=",0.019))


plot(mod$CCL20,mod$BMS_update_norm,mian="CCL20 with BMS sig.")
abline(lm(mod$BMS_update_norm~mod$CCL20),col="red")
# cor.test(mod$CCL20,mod$Brain_gene_norm,method="spearman")
legend("topright",legend=paste0("rho=",0.26," pvalue=",0.037))


plot(mod$CCL20,mod$Treg_norm,mian="CCL20 with BMS sig.")
abline(lm(mod$Treg_norm~mod$CCL20),col="red")
# cor.test(mod$CCL20,mod$Treg_norm,method="spearman")
legend("topright",legend=paste0("rho=",0.13," pvalue=",0.30))


dev.off()





#====================================================================================================================================================
# CCL20 in LUAD and LCBM 
# 2021-9-26

dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
tcga <- readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")

# housing-keeping 
grep("ACTB",rownames(dat),value=T)
grep("ACTB",rownames(tcga),value=T)

dat.f <- apply(dat,2,function(x){x/x[149]})
tcga.f <- apply(tcga,2,function(x){x/x[8176]})


pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/CCL20_LUAD_LCBM.pdf",useDingbats=F)
boxplot(dat.f["CCL20",],tcga.f["CCL20",],outline=F,names=c("LCBM","LUAD"),main="CCL20 (ACTB adjust)")
dev.off()




















#=============================================================================================================================================
# CCR6 in T cell subset expression percentage 
library(Seurat)
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS") # use all cell type 
dat <- subset(tmp,cells=which(tmp$celltype%in%c("CD4 naive","CD8 cytotoxic","CD8 exhausted/cytotoxic","Tfh","Treg")))
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$celltype)

features <- "CCR4"
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

res.dat.list <- percent_feature(dat,features,group="celltype")

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/Tcell_subset_CCR6.pdf",useDingbats=F)
barplot(sort(res.dat.list[[1]]),las=2)
dev.off()





#=======================================================================================================
# TFs activity in different group cells.
# 2021-9-24 
# plot in Fig4.R 
# 2021-12-4 
# use new cell type group to calculate 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
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

tf.dat.f$BMS.group <- dat$BMS.group
tmp <- aggregate(CEBPD~BMS.group,data=tf.dat.f,FUN=mean)

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/CEBPD_TFs.pdf",useDingbats=F)
barplot(c(0.11435307,0.08096245),names=c("BMS-high","BMS-low"),ylim=c(0,0.2),main="CEBPD activity")
dev.off()









##########################################################################################################################################################
# 2021-9-26
# do some verify for Treg enrichment 
# 0. check by GSE131907 data   # 2021-9-26 result show Treg percentage is not high in BM 
# 1. check in 123902
#=========================================================================================================================================================
# 0. GSE131907
library(Seurat)
tmp <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
dat <- subset(tmp,cells=which(tmp$Cell_type=="T lymphocytes"))









library(Seurat)
tmp <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")
DefaultAssay(tmp) <- "RNA"
FeaturePlot(tmp,features=c("CD3D","CD3E"),label=T)
VlnPlot(tmp,features=c("CD3D","CD3E"),pt.size=0)
tmp.dat <- subset(tmp,cells=which(tmp$seurat_clusters%in%c(0,1,3,5,8,9,12,17)))

# re-integration 
inte.list <- list()
samples <- unique(tmp.dat$sample)
for(i in 1:length(samples)){
    tmp <- subset(tmp.dat,cells=which(tmp.dat$sample==samples[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=50,k.score=50)
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

DefaultAssay(inte) <- "RNA"
FeaturePlot(inte,features=c("FOXP3","IL2RA"))
VlnPlot(inte,features=c("FOXP3","IL2RA"),pt.size=0)



# get Tregs 
FeaturePlot(tmp.dat,features=c("FOXP3","IL2RA"),label=T)





















#==============================================================================================================================================
# 2021-10-4
# GSE74639 analysis CTC with CCL20 
library(readxl)
tmp <- read_excel("/public/workspace/lily/metastasis/data/verify/GSE74639/GSE74639_readCounts.xls")
ann <- read_excel("/public/workspace/lily/metastasis/data/verify/GSE74639/GSE74639_annotation_for_readCounts.xls")
colnames(tmp)[1] <- "ID"
tmp.dat <- merge(tmp,ann[,c(1,2,4)],by="ID")
colnames(tmp.dat)[18] <- "Entrez.ID"
# tmp.dat$ID <- NULL

# counts to expreesion translation
anno <- read.table("~/REF/hg19_gene_ann.txt",sep="\t",header=T)
anno$length <- anno$Gene.end..bp. -anno$Gene.start..bp.
tmp.res <- merge(anno[,c(6,7)],tmp.dat,by.x="NCBI.gene..formerly.Entrezgene..ID",by.y="Entrez.ID")
tmp.res$NCBI.gene..formerly.Entrezgene..ID <- NULL
tmp.res$ID <- NULL
# transform into expression data 
################################################################################################################
# tmp.res[,2] <- NULL
# 
tmp <- t(apply(tmp.res[,-18],1,function(x){x/x[1]}))*10^3
res.f <- apply(tmp,2,function(x){x/sum(x)})*10^6
res.f <- data.frame(res.f)
res.f$gene_name <- tmp.res$symbol
res.data <- aggregate(.~gene_name,data=res.f,FUN=median)
rownames(res.data) <- res.data$gene_name
res.data[,c(1,2)] <- NULL
saveRDS(as.matrix(res.data),file="/public/workspace/lily/metastasis/data/verify/GSE74639/GSE74639_exp.RDS")

#=============================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE74639/GSE74639_exp.RDS")



############################################################################################################################
############################################################################################################################
############################################################################################################################
# 1. use GSE123902 DTC data 
# 2. use GSE123902 early parimary and metastasis 


# 1. GSE123902 DTC data 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro.RDS")
DefaultAssay(dat) <- "RNA"
dat$CCL20 <- ifelse(dat[["RNA"]]@data["CCL20",]>0,"Y","N")

tmp <- apply(table(dat$CCL20,dat$group),2,function(x){x/sum(x)})

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/CCL20_GSE123902_DTC_percent.pdf",useDingbats=F)
barplot(tmp[2,],main="GSE123902.DTC",ylab="CCL20 percentage")
dev.off()



# 2. GSE123902 early primary and advanced
# 2021-10-5
GSE123902 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")
GSE123902$group <- "Early"
GSE123902$group[which(GSE123902$sample=="LX675")] <- "Advanced"

GSE123902$CCL20 <- ifelse(GSE123902[["RNA"]]@data["CCL20",]>0,"Y","N")

tmp <- apply(table(GSE123902$CCL20,GSE123902$group),2,function(x){x/sum(x)})

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/CCL20_GSE123902_Primary_percent.pdf",useDingbats=F)
barplot(tmp[2,],main="GSE123902.Primary",ylab="CCL20 percentage")
dev.off()






















#============================================================================================================================================
# 2021-12-16
# use E-MTAB calculate CCL20 with BMS 

plot(mod$BMS_update_norm,as.numeric(dat.BM["CCL20",]),mian="CCL20 with BMS sig.",xlab="BMS score",ylab="Expression value of CCL20")
abline(lm(as.numeric(dat.BM["CCL20",])~mod$BMS_update_norm),col="red")
legend("topleft",legend=paste0("rho=",0.26," pvalue=",0.037),bty="n")
















#===========================================================================================================================================
# 2021-12-16
# check CCL20 expression in malignant cells and Myeloid cells
#===========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
dat.sub <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type=="maliganant"))

mye <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
mye.sub <- subset(mye,cells=names(which(mye$type_group=="LCBM"&mye$celltype%in%c("MDM","Monocyte","MG"))))

tmp.dat <- merge(dat.sub,mye.sub)
# re-integration 
inte.list <- list()
samples <- unique(tmp.dat$orig.ident)
for(i in 1:length(samples)){
    tmp <- subset(tmp.dat,cells=which(tmp.dat$orig.ident==samples[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}

integration.anchors <- FindIntegrationAnchors(inte.list)
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

########################################################################################
# check result 
inte$CCL20 <- ifelse(inte[["RNA"]]@data["CCL20",]>0,"Y","N")
FeaturePlot(inte,features="RNA_CCL20",order=T)
DimPlot(inte,group.by="type")

barplot(c(0.536752,0.4453288),names=c("malignant","Myeloid"),ylab="Number of CCL20+ cells",ylim=c(0,0.6))
















#=============================================================================================================================================
# 2021-12-17
# GSE59831 analysis to verify MDM expression CCL20 
#=============================================================================================================================================
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE59831/GSE59831_processed_data_FPKM.txt",sep="\t",header=T)

tmp.sample <- read.table("/public/workspace/lily/metastasis/data/verify/GSE59831/sampleinfo.txt",sep="\t")
tmp.sample$type <- sapply(strsplit(as.vector(tmp.sample$V3),": "),function(x){x[[2]]})

# use macrophage and monocyte 
gene <- "CCL20"
exp1 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("Tum1","Tum2","Tum3","Tum4","Tum5")])
exp2 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("WT1","WT2","WT3","WT4")])
exp3 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("Tum9","Tum10","Tum11")])
exp4 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("WT8","WT9","WT10")])
# exp5 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("Tum12","Tum13","Tum14")])
# exp6 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol==gene),c("Tum11","Tum12","Tum13")])


boxplot(exp1,exp2,
    names=c("Myeloid with tumor","Control Myeloid"),main="Myeloid CCL20"
)

boxplot(exp1,exp2,exp3,exp4,main="CCL20 expression",ylab="Expression",xlab="group",
    names=c("T.Mye.","C.Mye.","T.Epi.","C.Epi.")
)


# plot result 

tmp.sample$CCL20 <- as.numeric(tmp.dat[which(tmp.dat$human_gene_symbol=="CCL20"),-c(1,2)])

pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE59831_CCL20_myeloid.pdf",useDingbats=F)
boxplot(CCL20~V2,data=tmp.sample[c(1:9),],names=c("Myeloid with tumor","Control Myeloid"))
beeswarm::beeswarm(CCL20~V2,data=tmp.sample[c(1:9),],col ="black", pch = 19,main = 'beeswarm + bxplot',add = TRUE)
legend("topright",legend="P = 0.11")
dev.off()



pdf("/public/workspace/lily/Lung2Brain/Version5/T_cell/plot/GSE59831_CCL20_Epithelial.pdf",useDingbats=F)
boxplot(CCL20~V2,data=tmp.sample[c(16:21),],names=c("Tumor epithelial","Control Epithelial"))
beeswarm::beeswarm(CCL20~V2,data=tmp.sample[c(16:21),],col ="black", pch = 19,main = 'beeswarm + bxplot',add = TRUE)
legend("topright",legend="P = 0.1")
dev.off()

























#===========================================================================================
# do some trajectory for tumor 
# 2020-12-3
#===========================================================================================
library(Seurat)
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(Seurat)
library(monocle)

# first to find Variable Features to run monocle 
# 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$maliganant=="tumor"))

inte.list <- list() 
samplelist <- unique(sub.dat$orig.ident)
# samplelist <- samplelist[-7] #2020-12-15 change 
for(i in 1:length(samplelist)){
	tmp <- subset(sub.dat,cells=which(sub.dat$orig.ident==samplelist[i]))
    tmp_data <- tmp[["RNA"]]@data       #2020-12-15
    tmp.subset <- subset(tmp,cells=colnames(tmp_data)[which(tmp_data["MYBL2",]>0|tmp_data["CEBPB",]>0)])       #2020-12-15
	DefaultAssay(tmp.subset) <- "RNA"       #2020-12-15
	inte.list[i] <- tmp.subset        #2020-12-15 change 
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=20)
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
inte <- RunUMAP(inte,dims=1:10)

#==============================================================================================================
# transform data 
data <- as(as.matrix(inte@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = inte@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#====================================================================================
# get order genes 
#
#====================================================================================
DefaultAssay(inte) <- "integrated"
genes <-  rownames(inte)[1:1000]
#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- genes
cds <- setOrderingFilter(cds, ordering_genes)
# pseudotime 
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)

# plot_cell_trajectory(cds,color_by="type_group")+scale_colour_manual(values=c("#bff199","#70b29c"))+facet_wrap(~type_group, nrow = 3)

# # try to check TFs expression to show 
# cds$CEBPB <- cds@assayData$exprs["CEBPB",]
# plot_cell_trajectory(cds,color_by="CEBPB")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

# cds$MYBL2 <- cds@assayData$exprs["MYBL2",]
# plot_cell_trajectory(cds,color_by="MYBL2")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

# cds$MKI67 <- cds@assayData$exprs["MKI67",]
# plot_cell_trajectory(cds,color_by="MKI67")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

saveRDS(cds,file="/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory.RDS")
saveRDS(cds,file="/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS") # just CEBPB and MYBL2 expreesio cell 

#==============================================================================
# need to calculae BMS and TFs to show result 
# 2020-12-5 
#==============================================================================
# BMS calculate 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$maliganant=="tumor"))
#==============================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- data.frame(mod)
mod$type_group <- sub.dat$type_group
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_BMS.RDS")



# Pyscenic calculate 
# calculate in .shell script 
#=============================================================================================================================




# final show result 
# 2020-12-7
#=============================================================================================================================
# BMS mod score result plot 
library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$maliganant=="tumor"))
cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_BMS.RDS")
tf.dat <- read.delim2("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_auc_mtx.tsv",sep="\t",stringsAsFactors= F)
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))

#==============================================================================================================================
# plot result 
#==============================================
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/ALL_tumor_trajectory.pdf")
# 1. group by type_group
plot_cell_trajectory(cds,color_by="type_group")+scale_colour_manual(values=c("#bff199","#70b29c"))

# 2. group by TFs MYBL2 and CEBPB expression data 
cds$CEBPB_exp <- cds@assayData$exprs["CEBPB",]
plot_cell_trajectory(cds,color_by="CEBPB_exp")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

cds$MYBL2_exp <- cds@assayData$exprs["MYBL2",]
plot_cell_trajectory(cds,color_by="MYBL2_exp")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

cds$MKI67_exp <- cds@assayData$exprs["MKI67",]
plot_cell_trajectory(cds,color_by="MKI67_exp")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

# 3. group by TFs activaity
cds$MYBL2_ac <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
plot_cell_trajectory(cds,color_by="MYBL2_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

cds$CEBPB_ac <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
plot_cell_trajectory(cds,color_by="CEBPB_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))


# 4. group by MYBL2 and CEBPB expression bin
cds$MKI67_bin <- "no"
cds$MKI67_bin[which(cds@assayData$exprs["MKI67",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="MKI67_bin") + scale_colour_manual(values=c("#009f4d","#ff6c5f"))   

cds$CEBPB_bin <- "no"
cds$CEBPB_bin[which(cds@assayData$exprs["CEBPB",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="CEBPB_bin",) + scale_colour_manual(values=c("#009f4d","#ff6c5f"))    #+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.55,1))


cds$MYBL2_bin <- "no"
cds$MYBL2_bin[which(cds@assayData$exprs["MYBL2",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="MYBL2_bin")  + scale_colour_manual(values=c("#009f4d","#ff6c5f"))    # +scale_colour_gradientn(colors=c("steelblue","aliceblue","white",'orange'),values=c(0,0.4,0.55,1))


dev.off()

cds$MYBL2_ac <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
plot_cell_trajectory(cds,color_by="MYBL2_ac")+scale_colour_gradientn(colors=c('blue',"mediumpurple1","tan1",'orange'),values=c(0,0.45,1))

cds$CEBPB_ac <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
plot_cell_trajectory(cds,color_by="CEBPB_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))

scale_fill_gradientn(colours = c('steelblue',"white","white",'darkred'),
  values=c(0,MIN,MAX,1))

# pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_trajectory.pdf",useDingbats=F)
# plot_cell_trajectory(cds,color_by="CEBPB")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 
# plot_cell_trajectory(cds,color_by="MYBL2")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 
# plot_cell_trajectory(cds,color_by="MKI67")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 

# plot_cell_trajectory(cds,color_by="MYBL2_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))
# plot_cell_trajectory(cds,color_by="CEBPB_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))
# plot_cell_trajectory(cds,color_by="BMS")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.6,1))

# dev.off()





#========================================================================
# CEBPB and MYBL2 bin to show 
# 2020-12-17
# change colour 
#========================================================================
library(monocle)
library(Seurat)

cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/MYBL2_CEBPB_tumor_trajectory.pdf",useDingbats=F)
plot_cell_trajectory(cds,color_by="type_group")+scale_colour_manual(values=c("#205527","#0863b5"))
cds$type.TF <- "unknow"
cds$type.TF[which(cds@assayData$exprs["CEBPB",]>0&cds@assayData$exprs["MYBL2",]==0)] <- "CEBPB"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]==0)] <- "MYBL2"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]>0)] <- "Both"
plot_cell_trajectory(cds,color_by="type.TF")+scale_colour_manual(values=c("#ffbd0a","#46b7fd","#7ac142"))
dev.off()




#========================================================================================
# get BMS mod to show 
# 
#========================================================================================
# library(Seurat)
# library(monocle)
# cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory.RDS")
# mod <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_BMS.RDS")

# cds$BMS <- mod[,2]
# plot_cell_trajectory(cds,color_by="BMS")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.6,1))



#=============================================
# final to show 
# show 
#=============================================
# pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_trajectory.pdf",useDingbats=F)
# plot_cell_trajectory(cds,color_by="CEBPB")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 
# plot_cell_trajectory(cds,color_by="MYBL2")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 
# plot_cell_trajectory(cds,color_by="MKI67")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1)) # expression 

# plot_cell_trajectory(cds,color_by="MYBL2_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))
# plot_cell_trajectory(cds,color_by="CEBPB_ac")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))
# plot_cell_trajectory(cds,color_by="BMS")+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.6,1))

# dev.off()







#============================================================================================
# calculate LCBM cell BMS score 
# 2020-12-3
#============================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$maliganant=="tumor"&dat$type_group=="LCBM"))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("immuneEscape"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_tumor_MOD.RDS")





#==========================================================================================
# use Sangsung data to verify 
# 2020-12-3
#==========================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")
sub.dat <- subset(dat,cells=which(dat$Cell_subtype%in%c("Malignant cells","tS1","tS2","tS3")))
sub.dat.f <- subset(sub.dat,cells=which(sub.dat$type%in%c("tumor_Lung_advanced","tumor_Lung_early")))
saveRDS(sub.dat.f,file="/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_tumor.RDS")



# calculate BMS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat.f[["RNA"]]@data),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- data.frame(mod)


dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_tumor.RDS")
dat <- FindVariableFeatures(object = dat)
# scaling
all.genes <- rownames(x = dat)
tmp_dat <- ScaleData(object = dat, features = all.genes)





count_dat <- read.delim2("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt",sep="\t",header=T)
cellinfo <- read.delim2("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt",sep="\t",header=T)

rownames(count_dat) <- count_dat$Index
count_dat$Index <- NULL
rownames(cellinfo) <- cellinfo$Index
tmp <- CreateSeuratObject(
counts=count_dat,
project = "SeuratProject",
assay = "RNA",
min.cells = 0,
min.features = 0,
names.field = 1,
names.delim = "_",
meta.data = cellinfo
)

dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_tumor.RDS")
dat -> tmp_dat
tmp_dat = NormalizeData(object = tmp_dat)
tmp_dat <- FindVariableFeatures(object = tmp_dat)
# scaling
all.genes <- rownames(x = tmp_dat)
tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)


# do trajectory 
#==========================================================================================
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(Seurat)
library(monocle)

dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_tumor.RDS")
dat <- FindVariableFeatures(object = dat)

#====================================================================================
data <- as(as.matrix(dat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = dat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#====================================================================================
# get order genes 
#
#====================================================================================
# Find vvariable genes 


genes <-  VariableFeatures(dat)
#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- genes
cds <- setOrderingFilter(cds, ordering_genes)
# pseudotime 
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds,color_by="type")    +scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.6,1))


cds$MKI67 <- cds@assayData$exprs["MKI67",]
cds$MKI67 <- "no"
cds$MKI67[which(cds@assayData$exprs["MKI67",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="MKI67") + scale_colour_manual(values=c("#009f4d","#ff6c5f"))     #+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.45,1))


cds$CEBPB <- "no"
cds$CEBPB[which(cds@assayData$exprs["CEBPB",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="CEBPB",) + scale_colour_manual(values=c("#009f4d","#ff6c5f"))    #+scale_colour_gradientn(colors=c('blue','orange'),values=c(0,0.55,1))


cds$MYBL2 <- "no"
cds$MYBL2[which(cds@assayData$exprs["MYBL2",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="MYBL2")  + scale_colour_manual(values=c("#009f4d","#ff6c5f"))    # +scale_colour_gradientn(colors=c("steelblue","aliceblue","white",'orange'),values=c(0,0.4,0.55,1))







































#=================================================================================================================================================
# T cell analysis in Sangsung data 
# 1. Treg enrich in BM 
# 2. CD4 ANXA1 T cell 
#=================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")









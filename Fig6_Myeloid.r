
# 2021-6-7
# analysis Myeloid cells 
#======================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_myeloid.RDS")
inte.list <- list()
samplelist <- unique(dat$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_myeloid.RDS")





##########################################################################################################################################
# 2021-6-7
# define cell type 
#=========================================================================================================================================
library(Seurat)
library(SingleR)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_myeloid.RDS")
ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster<-SingleR(test=as.matrix(dat@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=dat$seurat_clusters,method="cluster")


# VlnPlot show result 
#========================================================================================================================================
# define MDM and MG and Monocyte 
# 2021-6-7
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig6/Data/inte9_myeloid.RDS")
dat$hbmarker <- NULL
dat$llymarker <- NULL
dat$celltype <- NULL
dat$putativeTumor3 <- NULL
dat$res.2 <- NULL
dat$Phase <- NULL
dat$G2M.Score <- NULL
dat$S.Score <- NULL
dat$percent.mito <- NULL
dat$nUMI <- NULL
dat$nGene <- NULL
dat$RNA_snn_res.2 <- NULL

dat$celltype <- "Myeloid"
dat$celltype[which(dat$seurat_clusters%in%c(14))] <- "Undefined"
dat$celltype[which(dat$seurat_clusters%in%c(1,7))] <- "Monocyte"
sub.dat <- subset(dat,cells=which(dat$celltype=="Myeloid"))
# calculate BMS score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/Myeloid_MDM_MG_mod.RDS")




#===========================================================================================================================
# define cell type 
mod$type <- "MDM"
mod$type[which(mod$BMDM_marker_norm<mod$MG_marker_norm)] <- "MG"
# table(mod$type)
all(names(dat$celltype[which(colnames(dat)%in%rownames(mod))])==rownames(mod))
dat$celltype[which(colnames(dat)%in%rownames(mod))] <- mod$type

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")











########################################################################################################################################
# plot result 
# 2021-6-8
# 1. TCGA purify
# 2. inte9 cell type percentage different 
# 3. DimPlot show MDM and MG 
# 4. LCBM and MDM 
# 5. CSOmap
# 6. Cellchat 
#=======================================================================================================================================
# 1. TCGA purify 
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/LUAD_absolute.txt",sep="\t",header=T)
tmp.dat$SampleID <- substr(gsub("-",".",tmp.dat$Sample.ID),1,15)
tmp.res <- aggregate(ABSOLUTE~SampleID,data=tmp.dat,FUN=median)
luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")

rownames(tmp.res) <- tmp.res$SampleID
res.f <- merge(tmp.res,luad_mod,by="row.names")


library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("~/Lung2Brain/Version5/Myeloid/TCGA_purity_BMS_cor.pdf",useDingbats=F)
p<-ggplot(res.f,aes(x=BMS_update_norm,y=ABSOLUTE)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Purity") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "spearman",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()






#=============================================================================================================================================
# 2021-11-22
# use LCBM data to do correlation with immune score 

dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GSE14108_estimate_score.gct",skip = 2,header = T)
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$purity <- as.numeric(as.vector(scores[4,-c(1,2)]))
mod$immune <- as.numeric(as.vector(scores[2,-c(1,2)]))


library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("~/Lung2Brain/Version5/Myeloid/GSE14108_purity_BMS_cor.pdf",useDingbats=F)
p<-ggplot(mod,aes(x=BMS_update_norm,y=immune)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("ESTIMATE Immune score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "spearman",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()























# 2. inte9 cell type 
# 
library(Seurat)
dat <- readRDS("~/Lung2Brain/Version5/Data/inte9_ntumor.RDS")
dat$type.g <- dat$type
dat$type.g[which(dat$type.g=="B_cell")] <- "lymphocyte"
dat$type.g[which(dat$type.g=="T_cell")] <- "lymphocyte"

tmp <- table(dat$type_group,dat$type.g)
tmp.f <- tmp[,-grep("malignant|unknow",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
##################################################################################
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,level=c("GBM","LCBM","LC"))
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/Cell_type_inte.pdf",useDingbats=F)
cols <- c('#377EB8','#910241','#B2DF8A','#F29403','#984EA3')
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()









# 3. DimPlot show MDM MG and Monocyte 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/Myeloid_DimPlot.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype",cols=c("#fbb034","#00a4e4","#89ba16","#caccd1"),reduction="tsne")
dev.off()

# and show cell type percentage 
tmp <- table(dat$type_group,dat$celltype)
tmp.f <- tmp[,-grep("malignant|Undefined",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,level=c("GBM","LCBM","LC"))
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/Myeloid_celltype.pdf",useDingbats=F)
cols <- c("#fbb034","#00a4e4","#89ba16","#caccd1")
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()





# 4. TCGA GBM and LCBM MDM and MG calculate 
#=================================================================================================================================
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("MG_marker","BMDM_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/GSE14108_LCBM_mod.RDS")

dat.tcga <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/genomicMatrix",sep="\t",header=T)
rownames(dat.tcga) <- dat.tcga$sample
dat.tcga$sample <- NULL
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.tcga),c("MG_marker","BMDM_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Myeloid/TCGA_GBM_mod.RDS")

#==========================================================================================
# plot result 
library(ggplot2)
library(reshape2)
tmp1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GSE14108_LCBM_mod.RDS")
dat.1 <- melt(tmp1[,c("MG_marker_norm","BMDM_marker_norm")])

# plot result2
tmp2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/TCGA_GBM_mod.RDS")
dat.2 <- melt(tmp2[,c("MG_marker_norm","BMDM_marker_norm")])

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GBM_LCBM_MG_MDM.pdf",useDingbats=F)
ggplot(dat.2, aes(x=value, fill=variable)) + geom_density()+theme_bw()+xlab("GBM")
ggplot(dat.1, aes(x=value, fill=variable)) + geom_density()+theme_bw()+xlab("LCBM")
dev.off()







# 5. CSOmap run tumor and Myeloid in LCBM and GBM 
#=================================================================================================================================
# 2021-6-8
# /public/workspace/lily/Lung2Brain/Version5/CSOmap/
# first is LCBM 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))
tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tumor$celltype <- "tumor"
tmp.res <- merge(sub.dat,tumor)
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
# maybe should remove Undefinded cells
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype%in%c("MDM","MG","Monocyte","tumor")))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")

# then is GBM 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="GBM"))
tmp1 <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
tumor.GBM <- subset(tmp1,cells=which(tmp1$type=="malignant"&tmp1$type_group=="GBM"))
tumor.GBM$celltype <- "GBM_tumor"
tmp.res <- merge(sub.dat,tumor.GBM)
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")

# last is LC 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Myeloid/inte9_myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="LC"))
tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LUAD_clust.RDS")
tumor$celltype <- "tumor"
tmp.res <- merge(sub.dat,tumor)
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")


# bytlib load r/4.0.4
.libPaths("/bioapps/Rlibs/4.0.4")
# load package
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo)

tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")
# function 1 
as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# function 2
run_CSOmap <- function(labelData, TPM, sampling=0, seed=315) {
    if (sampling) {
        TPM_ <- TPM
        labelData_ <- labelData

        set.seed(seed)
        index <- unlist(
            apply(as.data.frame(table(labelData_$labels) * sampling/ncol(TPM_)), 1, function(x) {
                sample(which(labelData_$labels==x[1]), x[2])
            })
        )

        TPM <- TPM_[,index]
        labelData <- labelData_[index,]
    }

    # Calculate optimized 3D coordinates
    affinityMat = getAffinityMat(TPM, LR, verbose = T)

    coords_res = runExactTSNE_R(
    X = affinityMat,
    no_dims = 3,
    max_iter = 1000,
    verbose = T
    )
    coords = coords_res$Y
    rownames(coords) <- colnames(TPM)
    colnames(coords) <- c('x', 'y', 'z')

    # Visualization(by 3D density)
    require(dplyr)
    # arrange data
    coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

    join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
    cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

    density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
    cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

    # p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")

    # Get significance
    signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
    contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)

    return(list(
        labelData=labelData,
        TPM=TPM,
        cellinfo_tbl=cellinfo_tbl,
        signif_results=signif_results,
        contribution_list=contribution_list
    ))

}

TPM <- as.matrix(tmp.res$RNA@counts)

labelData <- data.frame(
    cells  = colnames(tmp.res),
    labels = as.vector(tmp.res$celltype)
)

results = run_CSOmap(labelData, TPM, sampling=0)
saveRDS(results, file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_CSOmap.results.RDS")

saveRDS(results, file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_CSOmap.results.RDS")

saveRDS(results, file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_CSOmap.results.RDS")


# plot LCBM result 
#========================================================================================================================
library(plot3D)
result <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_CSOmap.results.RDS")
cellinfo_tbl <- result$cellinfo_tbl

stat = aggregate(cbind(x,y,z)~labels, cellinfo_tbl,mean)
rownames(stat) = stat[,1]
stat=stat[,-1]

cex=0.8
phi=30
theta=0
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='MDM'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#f2af00",
    xlim=range(cellinfo_tbl $x),ylim=range(cellinfo_tbl $y),
    zlim=range(cellinfo_tbl $z))   
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='tumor'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#ce1126",add=T)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='MG'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#5482ab",add=T)
tmp = cellinfo_tbl[which(cellinfo_tbl $labels=='Monocyte'),]
scatter3D(tmp$x, tmp$y, tmp$z,pch=16,cex=cex,
    phi= phi,theta= theta,col="#7ab800",add=T)

# plotly 
library(plotly)
set.seed(12345)
plot_ly(cellinfo_tbl, x = ~x, y = ~y, z = ~z, color = ~labels , text = ~labels,colors = c("#c9c9c9","#2db928","#fcd000","#d13814"),size=0.3)



#===================================================================================================================================
# plot Brplot to show distance 
result <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_CSOmap.results.RDS")
cellinfo_tbl <- result$cellinfo_tbl

stat = aggregate(cbind(x,y,z)~labels, cellinfo_tbl,mean)
rownames(stat) = stat[,1]
stat=stat[,-1]
dist(stat)
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_Myeloid_tumor_distance.pdf",useDingbats=F)
barplot(c(0.08978695,0.09794600,0.12931052),names=c("MDM","Monocyte","MG"),ylim=c(0,0.15))
dev.off()

#####################################################################
# GBM 
result <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_CSOmap.results.RDS")
cellinfo_tbl <- result$cellinfo_tbl

stat = aggregate(cbind(x,y,z)~labels, cellinfo_tbl,mean)
rownames(stat) = stat[,1]
stat=stat[,-1]
dist(stat)
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GBM_Myeloid_tumor_distance.pdf",useDingbats=F)
barplot(c(0.27407694,0.32412684,0.22030794),names=c("MDM","Monocyte","MG"),ylim=c(0,0.5))
dev.off()
############################################################################
# LUAD 
result <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_CSOmap.results.RDS")
cellinfo_tbl <- result$cellinfo_tbl

stat = aggregate(cbind(x,y,z)~labels, cellinfo_tbl,mean)
rownames(stat) = stat[,1]
stat=stat[,-1]
dist(stat)
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LUAD_Myeloid_tumor_distance.pdf",useDingbats=F)
barplot(c(0.4435861,0.5051479,0.6597093),names=c("MDM","Monocyte","MG"),ylim=c(0,0.8))
dev.off()




# plot together 
#=============================================================================
dat <- data.frame(percent=c(0.08978695,0.09794600,0.12931052,0.27407694,0.32412684,0.22030794,0.4435861,0.5051479,0.6597093),
    celltype=c(rep(c("MDM","Monocyte","MG"),3)),group=c(rep(c("LCBM","GBM","LUAD"),each=3)))


dat <- matrix(c(0.08978695,0.09794600,0.12931052,0.27407694,0.32412684,0.22030794,0.4435861,0.5051479,0.6597093),
    ncol=3,byrow=T)

colnames(dat) <- c("MDM","Monocyte","MG")
rownames(dat) <- c("LCBM","GBM","LUAD")

# radar chat is not good to use ; not suitable
# library(fmsb)
# radarchart( data.frame(dat)  , axistype=1 , 
#     #custom polygon
#     # pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
#     #custom the grid
#     cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,1,0.25), cglwd=0.8,
#     #custom labels
#     vlcex=0.8 
# )

library(ggplot2)
dat <- data.frame(distance=c(0.08978695,0.09794600,0.12931052,0.27407694,0.32412684,0.22030794,0.4435861,0.5051479,0.6597093),
    celltype=c(rep(c("MDM","Monocyte","MG"),3)),group=c(rep(c("LCBM","GBM","LUAD"),each=3)))

dat$celltype <- factor(dat$celltype,levels=c("MDM","Monocyte","MG"))

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/Myeloid_Tumor_distance.pdf",useDingbats=F)
ggplot(dat,aes(x=group,y=distance,fill=celltype,group=celltype)) + geom_bar(stat="identity",position="dodge")+ theme_classic()+
    scale_fill_manual(values=c("#fbb034","#89ba16","#00a4e4"))+
    geom_text(aes(label=round(distance,2),y = distance + 0.05) , position = position_dodge(0.9),size =3)
dev.off()













####################################################################################################################################
# run cell cell communication 
# maybe need to focus on BMDM ?
# 2021-6-8
library(Seurat)
library(CellChat)
# tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")
# tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")

#==================================================================================================================================
# Run CellChat 
data.input <- GetAssayData(tmp.res, assay = "RNA", slot = "data") # normalized data matrix
meta <- data.frame(tmp.res@meta.data) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
cellchat <- setIdent(cellchat, ident.use = "celltype") # set "labels" as default cell identity
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
saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/GBM.CellChat.RDS")

saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.CellChat.RDS")

saveRDS(cellchat,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/LUAD.CellChat.RDS")



##################################################################################################################################
# get cellchat 
# 2021-6-10
library(Seurat)
library(CellChat)

cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.CellChat.RDS")
groupSize <-  as.numeric(table(cellchat@idents))

# found LCBM MDM2tumor specific gene pair 
#===============================================================================================================================
# pdf("./tmp.pdf")
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# dev.off()
# MDM       MG Monocyte    tumor
# values <- c(1, 2, 3)
# names(values) <- c("p > 0.05", "0.01 < p < 0.05", "p < 0.01")
tmp <- netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
# df <- tmp[which(tmp$prob>0.05&tmp$pval>=3),]
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
#rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))
colnames(rs) <- c("MDM.tumor","MG.tumor","Monocyte.tumor")

# check MDM  specific interaction
a <- rs[which(rs$MG.tumor==0&rs$Monocyte.tumor==0),]
all(tmp[which(tmp$interaction_name_2%in%rownames(a)),"interaction_name_2"]==rownames(a))
a$pathway <- tmp[which(tmp$interaction_name_2%in%rownames(a)),"pathway_name"]

saveRDS(a,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/Pair/LCBM_MDM2tumor.specifc.RDS")




# Check GBM to have more speicfic 
cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/GBM.CellChat.RDS")
groupSize <-  as.numeric(table(cellchat@idents))
tmp <- netVisual_bubble(cellchat, sources.use = 2:4, targets.use = 1, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))
colnames(rs) <- c("MDM.tumor","MG.tumor","Monocyte.tumor")
saveRDS(rs,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/Pair/GBM_MDM2tumor.specifc.RDS")





# check LUAD to have more specific
cellchat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LUAD.CellChat.RDS")
groupSize <-  as.numeric(table(cellchat@idents))
tmp <- netVisual_bubble(cellchat, sources.use = 1:3, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","source.target","interaction_name_2")]
rs <- reshape2::dcast(tmp.res,source.target~interaction_name_2,value.var="prob")
rs[is.na(rs)] <- 0
rownames(rs) <- rs[,1]
rs <- rs[,-1]
rs <- data.frame(t(rs))
colnames(rs) <- c("MDM.tumor","MG.tumor","Monocyte.tumor")

a <- rs[which(rs$MG.tumor==0&rs$Monocyte.tumor==0),]
saveRDS(a,file="/public/workspace/lily/Lung2Brain/Version5/CellChat/Pair/LUAD_MDM2tumor.specifc.RDS")












#====================================================================================================================================
# 2021-6-17
# plot gene pair calculate by Cell_Cell_interaction.r 
#====================================================================================================================================
library(Seurat)
# merge cells and define cell group
dat1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
dat1$cell.group <- paste0(dat1$type_group,"_",dat1$celltype)
# all(names(dat1$cell.group[grep("LCBM_tumor",dat1$cell.group)])==rownames(dat1@meta.data[which(dat1$cell.group=="LCBM_tumor"),]))
dat1$cell.group[which(dat1$cell.group=="LCBM_tumor")] <- paste0("LCBM_tumor",dat1@meta.data[which(dat1$cell.group=="LCBM_tumor"),]$seurat_clusters)

dat2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")
dat2$cell.group <- paste0(dat2$type_group,"_",dat2$celltype)

dat3 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")
dat3$cell.group <- paste0(dat3$type_group,"_",dat3$celltype)

dat <- merge(dat1,y=c(dat2,dat3))


rece <- c("SDC1","SDC4","ITGA3","IL1R1","OSMR","EGFR")
ligand <- c("IL1B","OSM","THBS1","FN1","HBEGF")

# calculate 
dat@active.ident <- factor(dat$cell.group)

tumor <- AverageExpression(dat,assays="RNA",features=rece)$RNA
tumor <- tumor[,grep("tumor",colnames(tumor))]

mye <- AverageExpression(dat,assays="RNA",features=ligand)$RNA
mye <- mye[,grep("MDM",colnames(mye))]












# 2021-6-17
# GSE 137762 verify HBEGF 
#=======================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
Hbegf.res <- as.data.frame(dat["Ccl20",])
Hbegf.res$type <- "unknow"
colnames(Hbegf.res)[1] <- "exp"
Hbegf.res$type <- gsub("[0-9]$","",rownames(Hbegf.res))

#=======================================================================================================================
# just show Hbegf not show chemotherapy result 
#===========================================================
tmp.f <- Hbegf.res[grep("BMDM|Blood",rownames(Hbegf.res)),]
tmp.f <- tmp.f[1:14,]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM"))

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/HBEGF.BMDM.exp.pdf",useDingbats=F)
boxplot(exp~type,data=tmp.f,FUN=median)
beeswarm::beeswarm(exp~type,data=tmp.f,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()




# 2021-7-12
# MDM result change into IL1B and IL1R1 
#========================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
Il1b.res <- as.data.frame(dat["Ccl20",])
Il1b.res$type <- "unknow"
colnames(Il1b.res)[1] <- "exp"
Il1b.res$type <- gsub("[0-9]$","",rownames(Il1b.res))

#=======================================================================================================================
# just show Il1b not show chemotherapy result 
#===========================================================
tmp.f <- Il1b.res[grep("BMDM|Blood",rownames(Il1b.res)),]
tmp.f <- tmp.f[1:14,]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM"))

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/IL1B.BMDM.exp.pdf",useDingbats=F)
boxplot(exp~type,data=tmp.f,FUN=median)
beeswarm::beeswarm(exp~type,data=tmp.f,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()



#============================================================================================================================================
# 2021-12-4
# CEBPD verify 

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE137762/GSE137762_expr.RDS")
Cebpd.res <- as.data.frame(dat["Cebpd",])
Cebpd.res$type <- "unknow"
colnames(Cebpd.res)[1] <- "exp"
Cebpd.res$type <- gsub("[0-9]$","",rownames(Cebpd.res))

#=======================================================================================================================
# just show Cebpd not show chemotherapy result 
#===========================================================
tmp.f <- Cebpd.res[grep("BMDM|Blood",rownames(Cebpd.res)),]
tmp.f <- tmp.f[1:14,]
tmp.f$type <- factor(tmp.f$type,levels=c("Control_BloodMonocyte","small_BMDM","large_BMDM"))

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/Cebpd.BMDM.exp.pdf",useDingbats=F)
boxplot(exp~type,data=tmp.f,FUN=median,ylim=c(0,200))
beeswarm::beeswarm(exp~type,data=tmp.f,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()


















################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat)<- "RNA"
dat@active.ident <- factor(dat$group)

AverageExpression(dat,assays="RNA",features="IL1R1")$RNA

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/LCBM_subtumor_IL1R1.pdf",useDingbats=F)
barplot(c(3.336415,0.3584679,0.7695442),names=c("group09","group16","group47"),ylim=c(0,3.5))
dev.off()




# try to use GSE14108 to analysis
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)









#============================================================================================================================
# TCGA and GSE126548 to verify SDC4 
# 1. TCGA mod with SDC4 expression 
#============================================================================================================================
dat <- readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
mod$SDC4 <- as.numeric(as.vector(dat["SDC4",]))

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/TCGA.BMS.SDC4.pdf",useDingbats=F)
plot(mod$SDC4,mod$BMS_update_norm)
abline(lm(mod$BMS_update_norm~mod$SDC4),col="red")
legend("topright",legend=paste0("rho=",0.28))
dev.off()

# 2. GSE126548 to analysis show not significant result 


















#########################################################################################################################
# Analysis about Microgila
# 2021-6-21
#========================================================================================================================
library(Seurat)
library(CellChat)

cellchat1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LCBM.CellChat.RDS")
cellchat2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/GBM.CellChat.RDS")
# cellchat3 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CellChat/LUAD.CellChat.RDS")   # microglia of LUAD is not OK


# LCBM 
tmp <- netVisual_bubble(cellchat1, sources.use = 2, targets.use = 4, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "LCBM"
tmp1 <- tmp.res

# GBM 
tmp <- netVisual_bubble(cellchat2, sources.use = 3, targets.use = 1, remove.isolate = FALSE,return.data=T)$communication
tmp.res <- tmp[,c("prob","interaction_name_2")]
tmp.res$group <- "GBM"
tmp2 <- tmp.res

# merge data 
#===========================================================================
tmp.f <- merge(tmp1,tmp2,all.x=T,all.y=T,by="interaction_name_2")
colnames(tmp.f)[2:5] <- c("prob.LCBM","group.LCBM","prob.GBM","group.GBM")













#=======================================================================================================================================================
# 2021-9-23
# analysis about MDM and MG 
#=======================================================================================================================================================
# 1. use single cell data to analysis CCR6 and ITGB2 expression in MDM/MG 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$celltype)

AverageExpression(dat,assay="RNA",features="CCR6")


# 2. use bulk data to verify CCR6 with MDM/MG signature score 
# GSE14108
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


# E-MTAB
# use this because have CCL3 
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
















# 2021-12-3 
# use GSE161116 to do some verify
# however GSE161116 is not OK 
#==============================================================================================================================================
dat <- read.table("~/metastasis/data/verify/GSE161116/GSE161116_series_matrix.txt",sep="\t",header=T,comment.char="!")
rownames(dat) <- dat$ID_REF
dat$ID_REF <- NULL
info <- as.vector(read.table("~/metastasis/data/verify/GSE161116/sampleinfo.txt"))[-1]
names(info) <- NULL
ann <- data.frame(sample=colnames(dat),patient=sapply(sapply(info,function(x){strsplit(as.vector(x)," ")}),function(y){paste0(y[1],y[2])}),
    group=sapply(sapply(info,function(x){strsplit(as.vector(x)," ")}),function(y){paste0(y[3])})
)

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Treg"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)




#=============================================================================================================================================
# 2021-12-3
# use GSE14108 to calculate hallmark and IL1B ,IL1R1
dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat),modlist,'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- data.frame(mod)

tmp.res <- mod[,51:100]
tmp.res$IL1B <- as.numeric(dat["IL1B",])
tmp.res$IL1R1 <- as.numeric(dat["IL1R1",])

# calculate spearman
res.f <- data.frame(t(apply(tmp.res,2,function(x){
    c(
        cor.test(as.numeric(x),as.numeric(tmp.res[,"IL1B"]),method="spearman")$estimate,
        cor.test(as.numeric(x),as.numeric(tmp.res[,"IL1B"]),method="spearman")$p.value
    )
})))
colnames(res.f) <- c("Rho","pvalue")
res.final <- res.f[which(res.f$pvalue<0.05),]
res.final <- res.final[-c(22,23),]

p.data <- res.final[order(res.final$Rho,decreasing=T),]


#============================================================================================================================================
library(ggplot2)
p.data$Pathway <- gsub("^HALLMARK_|_norm$","",rownames(p.data))
p.data$Pathway <- factor(p.data$Pathway,levels=p.data$Pathway)
p.data$logP <-  -log2(p.data$pvalue)

pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GSE14108_IL1B_hallmark.pdf",useDingbats=F)
ggplot(p.data, aes(x = Pathway, y = Rho)) +
  geom_hline(yintercept = 0, color = "grey", size = 1) + # 添加y=0的辅助线
  geom_point(aes(color = Pathway,size=logP)) +         # 将点的size设置大一些比较好看
  geom_bar(aes(fill = Pathway), stat = "identity", width = 0.2) + # 注意将width宽度设小
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),      # 消除竖条的背景线
        axis.text.x = element_text(angle = 90),
        legend.position = "None",
        panel.border = element_blank(),
        # text = element_text(family = "STHeiti"), # Mac 电脑上绘图展现中文需要此行命令
        plot.title = element_text(hjust = 0.5)) 
dev.off()






# 2021-12-3 
# use GSE142620 to do verify
#=============================================================================================================================================
filelist <- dir()[-1]
tmp.res <- read.table(paste0("./",filelist[1]),header=F,sep="\t")
rownames(tmp.res) <- tmp.res$V1
colnames(tmp.res) <- c("ENSG","GSM4232997")
for(i in 2:length(filelist)){
    tmp <- read.table(paste0("./",filelist[i]),header=F,sep="\t")
    if(all(rownames(tmp.res)==tmp$V1)){
        tmp.res <- cbind(tmp.res,tmp[,2])
        colnames(tmp.res)[ncol(tmp.res)] <- strsplit(filelist[i],"_")[[1]][1]
    }
}

# now transform into TPM 
res.final <- TCount2TPM(tmp.res,"ENSG",36)

saveRDS(res.final,file="/public/workspace/lily/metastasis/data/verify/GSE142620/GSE142620_exp.RDS")

ann <- data.frame(samplename=colnames(res.final),treat=c(
    "Ctrl_6d_repA","IL1b_6d_repA","Ctrl_15d_repA","IL1b_15d_repA","Ctrl_21d_repA","IL1b_21d_repA","Ctrl_w6d_repA","IL1b_w6d_repA","Ctrl_w15d_repA","IL1b_w15d_repA","Ctrl_w30d_repA","IL1b_w30d_repA",
    "Ctrl_6d_repB","IL1b_6d_repB","Ctrl_15d_repB","IL1b_15d_repB","Ctrl_21d_repB","IL1b_21d_repB","Ctrl_w6d_repB","IL1b_w6d_repB","Ctrl_w15d_repB","IL1b_w15d_repB","Ctrl_w30d_repB","IL1b_w30d_repB",
    "Ctrl_6d_repC","IL1b_6d_repC","Ctrl_15d_repC","IL1b_15d_repC","Ctrl_21d_repC","IL1b_21d_repC","Ctrl_w6d_repC","IL1b_w6d_repC","Ctrl_w15d_repC","IL1b_w15d_repC","Ctrl_w30d_repC","IL1b_w30d_repC"
))

ann$group <- sapply(strsplit(as.vector(ann$treat),"_"),function(x){x[[1]]})
ann$time <- sapply(strsplit(as.vector(ann$treat),"_"),function(x){x[[2]]})
ann$rep <- sapply(strsplit(as.vector(ann$treat),"_"),function(x){x[[3]]})

saveRDS(ann,file="/public/workspace/lily/metastasis/data/verify/GSE142620/GSE142620_sampleinfo.RDS")

#=================================================================================================================================================
# now do analysis 
# 2021-12-3
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE142620/GSE142620_exp.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE142620/GSE142620_sampleinfo.RDS")

source('~/software/ssGSEA/ssgseaMOD.r')
# modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(as.matrix(dat),"BMS_update",'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

tmp.res <- mod[,2,drop=F]
tmp.res$group <- ann$group
tmp.res$time<- ann$time 
tmp.res$rep <- ann$rep 

#===========================================================================================================================================
tmp.res$CCL20 <- as.numeric(dat["CCL20",])
tmp.res$IL1R1 <- as.numeric(dat["IL1R1",])
tmp.res$IL1B <- as.numeric(dat["IL1B",])
tmp.res$CEBPD <- as.numeric(dat["CEBPD",])
tmp.res$class <- "before"
tmp.res$class[grep("w",tmp.res$time)] <- "after"

# 1. type : 
tmp.res$type <- paste0(tmp.res$group,"_",tmp.res$class)
tmp.res$type[grep("Ctrl",tmp.res$type)] <- "Ctrl"


# plot result 
#===========================================================================================================================================
pdf("/public/workspace/lily/Lung2Brain/Version5/Myeloid/GSE142620_IL1B_CCL20.pdf",useDingbats=F)
boxplot(CCL20~type,data=tmp.res,FUN=median)
beeswarm::beeswarm(CCL20~type,data=tmp.res,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)


boxplot(BMS_update_norm~type,data=tmp.res,FUN=median)
beeswarm::beeswarm(BMS_update_norm~type,data=tmp.res,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)


dev.off()




























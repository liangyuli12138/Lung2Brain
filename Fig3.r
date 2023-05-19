
#!/usr/bin/Rscript
# 2021-5-5
#=============================================================================================================
# 1. check TFs 
# 2. check High plascity 
# 3. trajecoty analysis



###############################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

# 0. update BMS score and High Plascity 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","HPSC_C5"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat$BMS <- NULL
dat$BMS_update <- mod[,3]
all(rownames(mod)==rownames(dat@meta.data))
dat$HPSC <- mod[,4]
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

# plot result
###############################
library(Seurat)
library(ggplot2)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/BMS_score_FeaturePlot.pdf",useDingbats=F)
DimPlot(dat,label=T)
FeaturePlot(dat,features="BMS_update") +  
    scale_colour_gradientn(colours =  c("#598c14","#b5c327","#edd812","#faae40","#ff3c41","red"),values = c(0,0.4,0.5,0.8,1.0))
dev.off()

# FeaturePlot(pbmc_small, c("LYZ", "MS4A1")) & 
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))




#================================================================================================================================
# Calculate Stmeness use different siganture 
# 2021-7-1
#================================================================================================================================
library(Seurat)
library(pheatmap)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
# calculate stemness signature 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("PNAS_stem","stem_cell"),"/public/workspace/lily/MOD_file/NatureMedicine/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_pur_stemness.RDS")


# result output 
library(pheatmap)
library(Seurat)
library(ggplot2)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_pur_stemness.RDS")
mod$BMS <- dat$BMS_update
dat.f <- mod[,c(3,4,7)]
dat.f$cytoTrace <- dat$cytoTrace

dat$PNAS_stem <- mod$PNAS_stem_norm
dat$cell_stem <- mod$stem_cell_norm

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/LCBM_pur_Stem.pdf",useDingbats=F)
FeaturePlot(dat,features=c("BMS_update"))+
scale_colour_gradientn(colours =  c("#598c14","#b5c327","#edd812","#faae40","#ff3c41","red"),values = c(0,0.4,0.5,0.8,1.0))

FeaturePlot(dat,features=c("cell_stem"))+
scale_colour_gradientn(colours =  c("#598c14","#b5c327","#edd812","#faae40","#ff3c41","red"),values = c(0,0.45,0.5,0.65,1.0))

pheatmap(dat.f,cluster_rows=F,cluster_cols=F,show_rownames=F,
    annotation_col=data.frame(row.names=colnames(dat.f),
        correlation=apply(dat.f,2,function(x){cor.test(x,dat.f$BMS,method="spearman")$estimate})
    ),
    color=colorRampPalette(c('steelblue','white','red'))(8000)
)

dev.off()





















# 1. Find high expression and activity TFs 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
dat$type.tumor <- "tumor.l"
dat$type.tumor[which(dat$seurat_clusters%in%c(4,7))] <- "tumor.h"  # identiy tumor subgroup
# 
tmp.data <- as.matrix(dat[['RNA']]@data)
tf.dat <- read.csv("/public/workspace/lily/Lung2Brain/Version5/Pyscenic/LCBM/step3.auc_mtx.csv")
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))

# identify gene expression
data.f <- tmp.data[which(rownames(tmp.data)%in%colnames(tf.dat)),]
index.h <- unname(which(dat$type.tumor=="tumor.h"))
tmp.res <- data.frame(t(apply(data.f,1,function(x){
    exp.h <- mean(x[index.h])
    exp.l <- mean(x[-index.h])
    log2fc <- log2(exp.h/exp.l)
    p <- wilcox.test(x[index.h],x[-index.h])$p.value
    c(exp.h,exp.l,log2fc,p)
})))
colnames(tmp.res) <- c("Exp.H","Exp.L","Log2FC","pvalue")
tmp.res$p.adj <- p.adjust(tmp.res$pvalue)

tmp.res <- tmp.res[order(tmp.res$Log2FC,decreasing=T),]
rs <- tmp.res[which(tmp.res$p.adj<0.05),]



# identify TF activity difference
# all(rownames(tf.dat)==colnames(dat))
# tmp.tf <- data.frame(t(apply(tf.dat,2,function(x){
#     exp.h <- mean(as.numeric(as.vector(x[index.h])))
#     exp.l <- mean(as.numeric(as.vector(x[-index.h])))
#     log2fc <- log2(exp.h/exp.l)
#     p <- wilcox.test(as.numeric(as.vector(x[index.h])),as.numeric(as.vector(x[-index.h])))$p.value
#     c(exp.h,exp.l,log2fc,p)
# })))
# colnames(tmp.tf) <- c("Exp.H","Exp.L","Log2FC","pvalue")
# tmp.tf$p.adj <- p.adjust(tmp.tf$pvalue)
# tmp.tf <- tmp.tf[order(tmp.tf$Log2FC,decreasing=T),]
# tf.rs <- tmp.tf[which(tmp.tf$p.adj<0.05),]


tf.dat$type.tumor <- dat$type.tumor
tmp <- aggregate(.~type.tumor,data=tf.dat,FUN=mean)
rownames(tmp) <- tmp$type.tumor
tmp$type.tumor <- NULL
tf.rs <- data.frame(t(tmp))
tf.rs$log2FC <- log2(tf.rs$tumor.h/tf.rs$tumor.l)
tf.rs <- tf.rs[order(tf.rs$log2FC,decreasing=T),]
# pheatmap::pheatmap(tmp)

# use Log2FC and log2FC to choose gene 
rs.f <- merge(rs,tf.rs,by="row.names")
rownames(rs.f) <- rs.f[,1]
colnames(rs.f)[1] <- "Gene"
rs.f$name <- rs.f$Gene
rs.f$name[-which(rs.f$Log2FC>0.5&rs.f$log2FC>0.5)] <- ""
# change Pvalue
rs.f$pvalue[which(rs.f$pvalue==0)] <- min(rs.f$pvalue[which(rs.f$pvalue>0)])
rs.f$log10Pvalue <- -log10(rs.f$pvalue)

library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/TFs.pdf",useDingbats=F)
ggplot(rs.f,aes(y=log10Pvalue,x=Log2FC,color=log2FC)) + geom_point(aes(size=2))+ geom_text(aes(label = name),nudge_x=0,nudge_y=0)+
    geom_vline(aes(xintercept=3.45)) + geom_vline(aes(xintercept=5.09)) +
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.8,0.8,1.0))+
    theme_bw()
dev.off()







# 2. plot E2F8 expression in LCBM clusters and in LC 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
# add a new meta data info 
dat$type.tumor <- "tumor.l"
dat$type.tumor[which(dat$seurat_clusters%in%c(4,7))] <- "tumor.h"  # identiy tumor subgroup
dat@active.ident <- factor(dat$type.tumor)
tmp <- AverageExpression(dat,assays="RNA",features="E2F8")



# LC <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LUAD_clust.RDS")
# LC@active.ident <- factor(LC$type_group)
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(LC[['RNA']]@data),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)
# median(mod[,2])
# AverageExpression(LC,assays="RNA",features="E2F8")


tmp.dat <- data.frame(Cluster=c("Tumor.H","Tumor.L"),stringsAsFactors=F,
    E2F8=as.numeric(unname(as.numeric(tmp$RNA))),
    BMS_update=as.numeric(aggregate(BMS_update~type.tumor,data=dat@meta.data,FUN=median)[,2]))
tmp.dat <- rbind(tmp.dat,c("LC",0.09190325,0.2979875))
tmp.dat$E2F8 <- as.numeric(tmp.dat$E2F8)
tmp.dat$BMS_update <- as.numeric(tmp.dat$BMS_update)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/E2F8_expression.pdf",useDingbats=F)
ggplot(tmp.dat,aes(x=Cluster,y=E2F8))+geom_bar(aes(fill=Cluster),stat="identity",position="dodge") + 
    theme_bw() + scale_fill_manual(values=c("blue","red","green"))
dev.off()



# 1. add a verify for E2F8 
# GSE14995
library(beeswarm)
GSE14995 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14995/GSE14995_dat.RDS")
tmp.dat <- data.frame(E2F8 = as.numeric(GSE14995["E2F8",]),group=c(rep("control",3),rep("invasion",3)))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/E2F8_GSE14995.pdf",useDingbats=F)
boxplot(E2F8~group,data=tmp.dat,main="E2F8 expression",names=c("Low invasion","High invasion"),ylim=c(0.9,1.2))
legend("topright",legend=(paste0("p=",wilcox.test(as.numeric(GSE14995["E2F8",1:3]),as.numeric(GSE14995["E2F8",4:6]))$p.value)))
beeswarm::beeswarm(E2F8~group,data=tmp.dat,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()














# 2. calculate high Plastic score 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
dat$type.tumor <- "tumor.l"
dat$type.tumor[which(dat$seurat_clusters%in%c(4,7))] <- "tumor.h"  # identiy tumor subgroup
aggregate(HPSC~type.tumor,data=dat@meta.data,FUN=median)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/LCBM_HPSC.pdf",useDingbats=F)
boxplot(HPSC~type.tumor,data=dat@meta.data,FUN=median,outline=F)
dev.off()










# 2021-5-14
# maybe do trajectory is not a good way to show result 
# 2021-8-12
# try to use trajectory to find specific genes 
# # 3. do a trajectory to show 
#=======================================================================================================================================
library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
# make a monocle object
DefaultAssay(dat) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(dat)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))

# find integration gene
inte.list <- list()
samplelist <- unique(dat$orig.ident)
for(i in 1:length(samplelist)){
    tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)

# run monocle 
DefaultAssay(inte) <- "integrated"
genes <- VariableFeatures(inte)[1:1000]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)

pdf("/public/workspace/lily/Lung2Brain/Version5/Trajectory/monocle/LCBM_subtumor_1000_monocle.pdf",useDingbats=F)
plot_cell_trajectory(tmp.dat,color_by="group")
plot_cell_trajectory(tmp.dat,color_by="BMS_update")
plot_cell_trajectory(tmp.dat,color_by="HPSC")
plot_cell_trajectory(tmp.dat,color_by="Lung.gene")
plot_cell_trajectory(tmp.dat,color_by="Brain.gene")
plot_cell_trajectory(tmp.dat,color_by="Pseudotime")
plot_cell_trajectory(tmp.dat,color_by="seurat_clusters") + facet_wrap(~seurat_clusters)
dev.off()

# saveRDS(tmp.dat,file="/public/workspace/lily/Lung2Brain/Version5/Trajectory/monocle/LCBM_tumor_monocle.RDS")

#===========================================================================================================================================
# try to find DEG for state 
# 2021-8-12
# use Seurat to find DEG of different State # because monocle DEG program result do not know how to check 
#===========================================================================================================================================
dat$State <- tmp.dat$State
dat@active.ident <- factor(dat$State)
gene <- FindAllMarkers(dat,min.pct=0,logfc.threshold=0)
gene <- gene[order(gene$avg_logFC,decreasing=T),]
























# # 4. veolocity to run 
# # bytlib load hdf5
# library(Seurat)
# library(loomR)

# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
# A05 <- subset(dat,cells=which(dat$orig.ident=="A20190305"))
# A12 <- subset(dat,cells=which(dat$orig.ident=="A20190312"))
# TB <- subset(dat,cells=which(dat$orig.ident=="T-Bsc1"))
# DefaultAssay(A05) <- "RNA"
# DefaultAssay(A12) <- "RNA"
# DefaultAssay(TB) <- "RNA"
# # Github say this to as.loom
# A05 <- FindVariableFeatures(A05)
# A12 <- FindVariableFeatures(A12)
# TB <- FindVariableFeatures(TB)

# A05$RNA_snn_res.2 <- NULL # debug for erro  https://github.com/mojaveazure/loomR/issues/40
# A12$RNA_snn_res.2 <- NULL
# TB$RNA_snn_res.2 <- NULL
# TB$RNA_snn_res.1 <- NULL

# # tansform into loom file 
# as.loom(A05,filename = "/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_tumor_veolcyto/A05.loom",verbose = T)
# as.loom(A12,filename = "/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_tumor_veolcyto/A12.loom",verbose = T)
# as.loom(TB,filename = "/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_tumor_veolcyto/TB.loom",verbose = T)









# 5. cytoTrace 
# R - 4.0.2

library(Seurat)
library(CytoTRACE)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

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

tmp.dat <- as_matrix(dat[['RNA']]@data)
# calculate 
results <- CytoTRACE(tmp.dat,ncores = 6,subsamplesize = 1000)
all(names(results$CytoTRACE) == colnames(dat))
dat$cytoTrace <- results$CytoTRACE
dat$type.tumor <- "tumor.l"
dat$type.tumor[which(dat$seurat_clusters%in%c(4,7))] <- "tumor.h" 
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
saveRDS(results,file="/public/workspace/lily/Lung2Brain/Version5/Data/cytotarce_LCBM_tumor.RDS")





# # GSE110495 data verfiy 
# # 2021-5-14
# # however do not know the exactlt meaning of BMIT LT BT BMIC 
# load("~/metastasis/data/verify/GSE110495/GSE110495.RData")
# load("~/metastasis/data/verify/GSE110495/GSE110495_anno.RData")

# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(gse110495),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)

# mod$type <- anno$stage.of.the.metastatic.cascade





# 6. plasticy and stemness result 
# plot Brain.gene signature and LC.gene signature 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tmp <- aggregate(cytoTrace~seurat_clusters,data=dat@meta.data,FUN=median)
tmp <- tmp[order(tmp$cytoTrace,decreasing=F),]
dat$seurat_clusters <- factor(dat$seurat_clusters,levels=as.vector(tmp$seurat_clusters))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Cytotrace_cluster.pdf",useDingbats=F)
boxplot(cytoTrace~seurat_clusters,data=dat@meta.data,FUN=median,outline=F)
dev.off()



# 7. plot Brain.gene and LC.gene signature 
# 2021-5-19 
# use circlize bar plot 
#=========================================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tmp <- aggregate(Lung.gene~seurat_clusters,data=dat@meta.data,FUN=median)
tmp <- tmp[-11,] # filter Cluster 10
colnames(tmp) <- c("Cluster","Lung.gene")
tmp$Cluster <- paste0("C",tmp$Cluster)
tmp <- tmp[order(tmp$Lung.gene,decreasing=T),]
tmp$Cluster <- factor(tmp$Cluster,levels=tmp$Cluster)

# Get the name and the y position of each label
label_data <- tmp
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (as.numeric(label_data$Cluster)-0.5) /number_of_bar     
# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
 
# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Tumor.Lung.gene.pdf",useDingbats=F)
library(ggplot2)
ggplot(data=tmp,aes(x=Cluster,y=Lung.gene)) + geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) + coord_polar(start = 0) + 
geom_text(data=label_data, aes(x=Cluster, y=Lung.gene+0.05, label=Cluster, hjust=hjust), color="black", fontface="bold", size=3) 
dev.off()













#  Brain.gene 
##################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
tmp <- aggregate(Brain.gene~seurat_clusters,data=dat@meta.data,FUN=median)
tmp <- tmp[-11,] # filter Cluster 10
colnames(tmp) <- c("Cluster","Brain.gene")
tmp$Cluster <- paste0("C",tmp$Cluster)
tmp <- tmp[order(tmp$Brain.gene,decreasing=T),]
tmp$Cluster <- factor(tmp$Cluster,levels=tmp$Cluster)

# Get the name and the y position of each label
label_data <- tmp
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (as.numeric(label_data$Cluster)-0.5) /number_of_bar     
# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
 
# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Tumor.Brain.gene.pdf",useDingbats=F)
library(ggplot2)
ggplot(data=tmp,aes(x=Cluster,y=Brain.gene)) + geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) + coord_polar(start = 0) + 
geom_text(data=label_data, aes(x=Cluster, y=Brain.gene+0.05, label=Cluster, hjust=hjust), color="black", fontface="bold", size=3) 
dev.off()





# maybe add a feature plot is more good
# 2021-5-19
#=======================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Lung.Brain.FeaturePlot.pdf",useDingbats=F)
FeaturePlot(dat,features="Lung.gene",min.cutoff=median(dat$Lung.gene),order=T)
FeaturePlot(dat,features="Brain.gene",min.cutoff=median(dat$Brain.gene),order=T)
dev.off()





# 2021-5-19
# maybe stemness need to calculate ?

# 2021-5-19
# use Lung to Brain double lesion to verify
#========================================================================================================
# D0927
library(Seurat) 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
# calculate signature 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat$BMS.update <- mod[,4]
dat$Brain.gene <- mod[,5]
dat$Lung.gene <- mod[,6]
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")

library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Fig3_D0927_Tumor.pdf",useDingbats=F)
DimPlot(dat)
FeaturePlot(dat,features="BMS.update",order=T)+
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.7,0.7,1.0))  
FeaturePlot(dat,features="Lung.gene",order=T) + 
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.6,0.6,1.0))  
FeaturePlot(dat,features="Brain.gene",order=T) +
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.45,0.45,1.0)) 
# DimPlot(dat,reduction="tsne")
# FeaturePlot(dat,features="BMS.update",min.cutoff=median(dat$BMS.update),reduction="tsne",order=T)
# FeaturePlot(dat,features="Lung.gene",min.cutoff=median(dat$Lung.gene),reduction="tsne",order=T)
# FeaturePlot(dat,features="Brain.gene",min.cutoff=median(dat$Brain.gene),reduction="tsne",order=T)
dev.off()




# E0927
library(Seurat) 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
# calculate signature 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat$BMS.update <- mod[,4]
dat$Brain.gene <- mod[,5]
dat$Lung.gene <- mod[,6]
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
# then recluster 
recluster <- function(tmp_dat){

# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=1)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)

	return(tmp_dat)
}



library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Fig3_E0927_Tumor.pdf",useDingbats=F)
DimPlot(dat)
FeaturePlot(dat,features="BMS.update",order=T)+
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.6,0.6,1.0))  
FeaturePlot(dat,features="Lung.gene",order=T) + 
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.6,0.6,1.0))  
FeaturePlot(dat,features="Brain.gene",order=T) +
    scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.55,0.55,1.0)) 
# DimPlot(dat,reduction="tsne")
# FeaturePlot(dat,features="BMS.update",min.cutoff=median(dat$BMS.update),reduction="tsne",order=T)
# FeaturePlot(dat,features="Lung.gene",min.cutoff=median(dat$Lung.gene),reduction="tsne",order=T)
# FeaturePlot(dat,features="Brain.gene",min.cutoff=median(dat$Brain.gene),reduction="tsne",order=T)
dev.off()







# use time data to veirfy if E2F8 is increas
# GSE136935
# not oK

library(reshape2)
library(dplyr)
library(GEOquery)
library(idmap1)
gset<-getGEO('GSE136935',destdir='/public/workspace/wangzy/work/LungMet/Data/')
gset=gset[[1]] 
Edata=exprs(gset)
map = idmap1::getIDs('GPL13497')
rownames(map)=map[,1]
Edata = merge(map,Edata,by='row.names')
Edata = as.data.frame(Edata[,-c(1:2,4)]%>%group_by(symbol)%>%summarise_all(mean))
Edata = Edata[!is.na(Edata[,1]),]
rownames(Edata)=Edata[,1]
Edata=Edata[,-1]
pd=pData(gset)
pd = pd[pd$`treatmnet:ch1`=='PBS',]
pd = pd[pd$`cell line:ch1`=='A549',]

####################################################################################################################
info <- pd[,c(33,34,35)]
colnames(info) <- c("cell_line","treattime","treatment")
all(rownames(info)==colnames(Edata))
dat <- Edata[,which(info$treatment=="PBS"&info$cell_line=="A549")]

# calculate BMS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


info.f <- info[which(info$treatment=="PBS"&info$cell_line=="A549"),]
tmp.res <- data.frame(E2F8=as.numeric(as.vector(dat["E2F8",])),time=as.vector(info.f$treattime))







# GSE58355
# plot E2F8 expression and BMS score 
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE58355/GSE58355_series_matrix.txt",comment.char="!",header=T,sep="\t")
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GSE58355/gene_ann.txt",sep="\t",header=T)
tmp.res <- merge(tmp.dat,ann,by.x="ID_REF",by.y="ID")

tmp.res.f <- tmp.res[-grep("\\|",tmp.res$Gene_symbol),]
tmp.res.f$ID_REF <- NULL
res.f <- aggregate(.~Gene_symbol,data=tmp.res.f,FUN=median)
rs <- res.f[-c(1,2),]

saveRDS(rs,file="/public/workspace/lily/metastasis/data/verify/GSE58355/GSE58355_exp.RDS")

#======================================================================================================================
# caculate BMS score 
# 2021-6-3

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE58355/GSE58355_exp.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[,c(1:6)]),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

library(beeswarm)
tmp.dat <- data.frame(BMS=mod[,2],group=c(rep("day2",3),rep("day25",3)))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/GSE58355_primary_BMS.pdf",useDingbats=F)
boxplot(BMS~group,data=tmp.dat,main="BMS score in primary tumor ",names=c("day2","day25"),ylim=c(0,1))
legend("topright",legend=(paste0("p=",wilcox.test(BMS~group,data=tmp.dat)$p.value)))
beeswarm::beeswarm(BMS~group,data=tmp.dat,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()




#=================================================================================================================
# calculate E2F8 
# 2021-6-3

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE58355/GSE58355_exp.RDS")

library(beeswarm)
tmp.dat <- data.frame(E2F8=as.numeric(as.vector(dat["E2F8",13:18])),group=c(rep("day10",3),rep("day25",3)))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/GSE58355_metastasis_E2F8.pdf",useDingbats=F)
boxplot(E2F8~group,data=tmp.dat,main="E2F8 in metastasis tumor ",names=c("day10","day25"))
legend("topright",legend=(paste0("p=",wilcox.test(E2F8~group,data=tmp.dat)$p.value)))
beeswarm::beeswarm(E2F8~group,data=tmp.dat,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)
dev.off()






#=================================================================================================================
# calculate E2F8 in BM- and BM+
# GSE




















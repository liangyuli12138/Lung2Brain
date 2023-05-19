
# this program is used to do tumor trajecoty 2021-5-13
# 1. do PATHE 
bytlib load phate
library(Seurat)

#.libPaths(c('/public/workspace/liangyuan/pipelines/bitapps_phate/ly_5873c054-b79d-478e-9ae5-bc6da1635829_RV3.6.0/','/bioapps/Rlibs/3.6.0/',"/public/workspace/lily/R/x86_64-pc-linux-gnu-library/3.6.0"))
#install.packages("phateR")

library(Seurat)
library(phateR)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
# filter some genes
mat<-as.matrix(dat@assays$RNA@counts)
keep_rows <- rowSums(mat>0) > 10
mat2<-mat[keep_rows,]
mat2 <- library.size.normalize(mat2)
mat2 <- sqrt(mat2)
pe<-phate(t(mat2),ndim=2,knn=20,decay=1000,gamma=1,t="auto")
dat[["phate"]] <- CreateDimReducObject(embeddings = pe$embedding, key = "phate_", assay = DefaultAssay(dat))
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Trajectory/PHATE/LCBM_tumor_phate.RDS")
# tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Trajectory/monocle/LCBM_tumor_monocle.RDS")




# 2. do DPT 

library(destiny)
library(Seurat)
library(Biobase)

knn=NULL
n_eigs=NULL
density_norm=NULL
distance_method=NULL
tips_method=NULL
seed_num=NULL

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
as.matrix(dat@assays$RNA@data)->exp
anno=dat@meta.data
es <- Biobase::ExpressionSet(exp, phenoData=Biobase::AnnotatedDataFrame(anno))

dmap <- DiffusionMap(es,sigma = "local",verbose=T,n_eigs=20,
                    density_norm = TRUE,distance = "euclidean")

# add a meta info 
dat$type.tumor <- "tumor.l"
dat$type.tumor[which(dat$seurat_clusters%in%c(4,7))] <- "tumor.h"
dpt <- DPT(dmap)





# 3. TSCAN 
# use R-4.0.2
library(TSCAN)
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
seuratdf<-as.matrix(dat@assays$RNA@counts)
procdata <- preprocess(seuratdf,cvcutoff = 0)
lpsmclust <- exprmclust(procdata)
mclustobj <- lpsmclust
saveRDS(mclustobj,file="/public/workspace/lily/Lung2Brain/Version5/Trajectory/TSCAN/LCBM_tumor_TSCAN.RDS")
# make a dataframe to plot result 
# lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid,sample_name = names(mclustobj$clusterid))
# S_matrix <- mclustobj$pcareduceres
# pca_space_df <- data.frame(S_matrix[, c(1, 2)])
# colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
# pca_space_df$sample_name <- row.names(pca_space_df)
# edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", by.y = "sample_name") # plot data 
# rownames(edge_df) <- edge_df$sample_name
# edge_df <- edge_df[colnames(dat),]
# all(rownames(edge_df)==colnames(dat))

# 2021-5-13 maybe just use pseudotime to plot 











































# 3.1 use double lung cacner brain metastasis sample 
library(Seurat) 
library(monocle)
dat.d <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.d[['RNA']]@data),c("BMS_update","HPSC_C5","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat.d$BMS_update <- mod[,5]
dat.d$HPSC <- mod[,6]
dat.d$Brain_gene <- mod[,7]
dat.d$Lung_gene <- mod[,8]

DefaultAssay(dat.d) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(dat.d)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))
dat.d <- FindVariableFeatures(dat.d)
genes <- VariableFeatures(dat.d)[1:500]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)
# plot_cell_trajectory(tmp.dat,color_by="State")

pdf("tmp_d.pdf")
plot_cell_trajectory(tmp.dat,color_by="BMS_update")
plot_cell_trajectory(tmp.dat,color_by="HPSC")
plot_cell_trajectory(tmp.dat,color_by="Lung_gene")
plot_cell_trajectory(tmp.dat,color_by="Brain_gene")
plot_cell_trajectory(tmp.dat,color_by="State")
dev.off()

aggregate(Lung_gene~State,data=pData(tmp.dat),FUN=median)


# another sample 
library(Seurat) 
library(monocle)
dat.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.e[['RNA']]@data),c("BMS_update","HPSC_C5","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat.e$BMS_update <- mod[,5]
dat.e$HPSC <- mod[,6]
dat.e$Brain_gene <- mod[,7]
dat.e$Lung_gene <- mod[,8]

DefaultAssay(dat.e) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(dat.e)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))
dat.d <- FindVariableFeatures(dat.e)
genes <- VariableFeatures(dat.e)[1:1500]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)
# plot_cell_trajectory(tmp.dat,color_by="State")

pdf("tmp_e.pdf")
plot_cell_trajectory(tmp.dat,color_by="BMS_update")
plot_cell_trajectory(tmp.dat,color_by="HPSC")
plot_cell_trajectory(tmp.dat,color_by="Lung_gene")
plot_cell_trajectory(tmp.dat,color_by="Brain_gene")
plot_cell_trajectory(tmp.dat,color_by="State")
dev.off()

aggregate(Lung_gene~State,data=pData(tmp.dat),FUN=median)


# use two sample integration result 
library(Seurat) 
library(monocle)
dat.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
dat.d <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")

integration.anchors <- FindIntegrationAnchors(object.list = c(dat.e,dat.d))
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
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.DE.RDS")

# run trajectory 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.DE.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","HPSC_C5","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat$BMS_update <- mod[,5]
dat$HPSC <- mod[,6]
dat$Brain_gene <- mod[,7]
dat$Lung_gene <- mod[,8]

DefaultAssay(dat) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(dat)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))
DefaultAssay(dat) <- "integrated"
genes <- VariableFeatures(dat)[1:1600]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)


pdf("tmp_DE.pdf")
plot_cell_trajectory(tmp.dat,color_by="BMS_update")
plot_cell_trajectory(tmp.dat,color_by="HPSC")
plot_cell_trajectory(tmp.dat,color_by="Lung_gene")
plot_cell_trajectory(tmp.dat,color_by="Brain_gene")
plot_cell_trajectory(tmp.dat,color_by="State")
dev.off()

aggregate(Lung_gene~State,data=pData(tmp.dat),FUN=median)




















































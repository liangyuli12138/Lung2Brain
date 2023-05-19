
# 2021-8-6
# do trajectory analysis 
#=========================================================================================================
# do trajectory for LCBM sub tumor 
# re-integration and do analysis
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")

# re-integration 
DefaultAssay(dat) <- "RNA"
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
##Scaling the integrateda
all.genes <- rownames(inte)
inte <- ScaleData(inte, features = all.genes)
#PCA
inte <- RunPCA(inte)
#cluster
inte <- FindNeighbors(inte)
inte <- FindClusters(inte,resolution=0.5)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)

#=====================================================================================================
# run monocle
library(monocle)
DefaultAssay(inte) <- "RNA"
tmp.dat <- Seurat::as.CellDataSet(inte)
tmp.dat <- estimateSizeFactors(tmp.dat)
tmp.dat <- detectGenes(tmp.dat, min_expr = 1)
# fData(dat)$use_for_ordering <-fData(dat)$num_cells_expressed > 0.05 * ncol(dat)
expressed_genes <- row.names(subset(fData(tmp.dat),num_cells_expressed >= 10))

DefaultAssay(inte) <- "integrated"
genes <- VariableFeatures(inte)[1:800]
ordering_genes <- genes
tmp.dat <- setOrderingFilter(tmp.dat, ordering_genes = ordering_genes)
tmp.dat <- reduceDimension(tmp.dat, method = 'DDRTree')
tmp.dat <- orderCells(tmp.dat)







































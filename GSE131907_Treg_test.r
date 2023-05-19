
#======================================================================================================
# 2020-11-30
# use samgsung data to run monocle 
#======================================================================================================
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(Seurat)
library(monocle)

dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")

#=======================================================================================================
tmp.dat <- subset(dat,cells=which(dat$Cell_subtype%in%c("Naive CD4+ T","Treg","CD4+ Th")))
tmp.dat.f <- subset(tmp.dat,cells=which(tmp.dat$type%in%c("Brian_metastasis","tumor_Lung_early")))


tmp_dat <- FindVariableFeatures(object = tmp.dat.f,nfeatures= 3000)
# scaling
all.genes <- rownames(x = tmp_dat)
tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)


#======================================================================================================
data <- as(as.matrix(tmp_dat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = tmp_dat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# make a object
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- VariableFeatures(tmp_dat)
cds <- setOrderingFilter(cds, ordering_genes)
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)


plot_cell_trajectory(cds,color_by="Cell_subtype")+scale_colour_manual(values=c("#bff199","#bff199","#70b29c"))+facet_wrap(~type, nrow = 3)

plot_cell_trajectory(cds,color_by="Cell_subtype")+scale_colour_manual(values=c("#bff199","#bff199","#70b29c"))+facet_wrap(~State, nrow = 3)


treg <- pData(cds)[which(cds$State%in%c(4,5)&cds$Cell_subtype=="Treg"),c("Cell_subtype","State")]

# treg.dat <- subset(tmp_dat,cells=which(colnames(tmp_dat)%in%rownames(treg)&tmp_dat$type%in%c("Brian_metastasis","tumor_Lung_early")))
treg.dat <- subset(tmp_dat,cells=which(colnames(tmp_dat)%in%rownames(treg)&tmp_dat$type%in%c("Brian_metastasis")))
treg.dat <- subset(tmp_dat,cells=which(colnames(tmp_dat)%in%rownames(treg)&tmp_dat$type%in%c("tumor_Lung_early")))
treg.dat$State <- treg$State[which(rownames(treg)%in%colnames(treg.dat))]


treg.dat@active.ident <- as.factor(treg.dat$State)
DefaultAssay(treg.dat) <- "RNA"
geneset_tmp <- FindMarkers(treg.dat,ident.1=4,ident.2=5,assay="RNA",logfc.threshold = 0.1)































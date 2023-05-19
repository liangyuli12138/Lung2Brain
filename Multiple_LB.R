

# 2021-4-6
# run multiple Lung cancer sample in Lung2Brain 
# run in 202.195.187.3
# D0927 and E0927
# 季桂枝 省人医
#====================================================================================================
library(Seurat)
# D0927
tmp <- Read10X(data.dir = "/public/workspace/lily/Mutiple_LB/D0927/D0927/outs/filtered_feature_bc_matrix/")
tmp.dat <- CreateSeuratObject(counts = tmp,  project = "D0927",min.cells = 3, min.features = 200)

prepare <- function(tmp_dat){
    tmp_dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp_dat, pattern = "^MT-")
    tmp_dat = subset(x=tmp_dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
# seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
	return(tmp_dat)
}
tmp.dat <- prepare(tmp.dat)
dat <- tmp.dat
# use genemarker 
    pdf(paste0("/public/workspace/lily/Mutiple_LB/D0927_violin.pdf"),width=12,height=9)
    DefaultAssay(dat) <- "RNA"
    VlnPlot(dat,features=c("CD3D","CD3E","CD2","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # T cells
    VlnPlot(dat,features=c("CD19","MS4A1","CD79A","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # B cells
    VlnPlot(dat,features=c("CD68","FCGR3A","LYZ","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Myeloid 
    VlnPlot(dat,features=c("NCAM1","NKG7","CD3D","KLRD1"),assay="RNA",group.by="seurat_clusters",pt.size=0) # NK cells
    VlnPlot(dat,features=c("COL1A1","COL1A2","THY1","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Fibroblast
    VlnPlot(dat,features=c("PTPRC","CLDN5","FLT1","RAMP2"),assay="RNA",group.by="seurat_clusters",pt.size=0) #Endothelials
    VlnPlot(dat,features=c("EPCAM","KRT19","CDH1","EGFR"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Epithelials
    VlnPlot(dat,features=c("TPSAB1","TPSB2","MS4A2","KIT"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Mast cells
    VlnPlot(dat,features=c("CLEC10A","CD1C", "CLEC4C", "IL3RA"),assay="RNA",group.by="seurat_clusters",pt.size=0) # pDC 
    dev.off()

    pdf(paste0("/public/workspace/lily/Mutiple_LB/D0927_feature.pdf"),width=10,height=10)
    DimPlot(dat,label=T)
    FeaturePlot(dat,features=c("CD3D","CD3E","CD2","PTPRC"),label=T) # T cells
    FeaturePlot(dat,features=c("CD19","MS4A1","CD79A","PTPRC"),label=T) # B cells
    FeaturePlot(dat,features=c("CD68","FCGR3A","LYZ","PTPRC"),label=T) # Myeloid 
    FeaturePlot(dat,features=c("NCAM1","NKG7","CD3D","KLRD1"),label=T) # NK cells
    FeaturePlot(dat,features=c("COL1A1","COL1A2","THY1","PTPRC"),label=T) # Fibroblast
    FeaturePlot(dat,features=c("PTPRC","CLDN5","FLT1","RAMP2"),label=T) #Endothelials
    FeaturePlot(dat,features=c("EPCAM","KRT19","CDH1","EGFR"),label=T) # Epithelials
    FeaturePlot(dat,features=c("TPSAB1","TPSB2","MS4A2","KIT"),label=T) # Mast cells
    FeaturePlot(dat,features=c("CLEC10A","CD1C","CLEC4C", "IL3RA"),label=T) # pDC 
    dev.off()

dat$infercnv.type <- "tmp"
dat$infercnv.type[which(dat$seurat_clusters%in%c(2,15))] <- "Tcell"
dat$infercnv.type[which(dat$seurat_clusters%in%c(0,1,4,8,11))] <- "Tumor"
saveRDS(dat,file="/public/workspace/lily/Mutiple_LB/D0927.RDS")


###############################################################################################################################
# run infercnv 
library(infercnv)
library(Seurat)
D0927 <- readRDS("/public/workspace/lily/Mutiple_LB/D0927.RDS")
sub.dat <- subset(D0927,cells=which(D0927$infercnv.type%in%c("Tcell","Tumor")))


setwd('/public/workspace/lily/Mutiple_LB/infercnv/')
dir.create('./D0927')
write.table(sub.dat$infercnv.type,file="D0927_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub.dat@assays$RNA@counts)
write.table(count,file="D0927_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/public/workspace/lily/Mutiple_LB/infercnv/D0927_count_exp.txt",
         annotations_file="/public/workspace/lily/Mutiple_LB/infercnv/D0927_cell_info.txt",
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names="Tcell")


setwd('/public/workspace/lily/Mutiple_LB/infercnv/D0927')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=10, #big
                             no_plot=F ,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )






#===========================================================================================================================
# E0927
tmp <- Read10X(data.dir = "/public/workspace/lily/Mutiple_LB/E0927/E0927/outs/filtered_feature_bc_matrix/")
tmp.dat <- CreateSeuratObject(counts = tmp,  project = "E0927",min.cells = 3, min.features = 200)

prepare <- function(tmp_dat){
    tmp_dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp_dat, pattern = "^MT-")
    tmp_dat = subset(x=tmp_dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
# seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
	return(tmp_dat)
}
tmp.dat <- prepare(tmp.dat)
dat <- tmp.dat
# use genemarker 
    pdf(paste0("/public/workspace/lily/Mutiple_LB/E0927_violin.pdf"),width=12,height=9)
    DefaultAssay(dat) <- "RNA"
    VlnPlot(dat,features=c("CD3D","CD3E","CD2","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # T cells
    VlnPlot(dat,features=c("CD19","MS4A1","CD79A","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # B cells
    VlnPlot(dat,features=c("CD68","FCGR3A","LYZ","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Myeloid 
    VlnPlot(dat,features=c("NCAM1","NKG7","CD3D","KLRD1"),assay="RNA",group.by="seurat_clusters",pt.size=0) # NK cells
    VlnPlot(dat,features=c("COL1A1","COL1A2","THY1","PTPRC"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Fibroblast
    VlnPlot(dat,features=c("PTPRC","CLDN5","FLT1","RAMP2"),assay="RNA",group.by="seurat_clusters",pt.size=0) #Endothelials
    VlnPlot(dat,features=c("EPCAM","KRT19","CDH1","EGFR"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Epithelials
    VlnPlot(dat,features=c("TPSAB1","TPSB2","MS4A2","KIT"),assay="RNA",group.by="seurat_clusters",pt.size=0) # Mast cells
    VlnPlot(dat,features=c("CLEC10A","CD1C", "CLEC4C", "IL3RA"),assay="RNA",group.by="seurat_clusters",pt.size=0) # pDC 
    dev.off()

    pdf(paste0("/public/workspace/lily/Mutiple_LB/E0927_feature.pdf"),width=10,height=10)
    DimPlot(dat,label=T)
    FeaturePlot(dat,features=c("CD3D","CD3E","CD2","PTPRC"),label=T) # T cells
    FeaturePlot(dat,features=c("CD19","MS4A1","CD79A","PTPRC"),label=T) # B cells
    FeaturePlot(dat,features=c("CD68","FCGR3A","LYZ","PTPRC"),label=T) # Myeloid 
    FeaturePlot(dat,features=c("NCAM1","NKG7","CD3D","KLRD1"),label=T) # NK cells
    FeaturePlot(dat,features=c("COL1A1","COL1A2","THY1","PTPRC"),label=T) # Fibroblast
    FeaturePlot(dat,features=c("PTPRC","CLDN5","FLT1","RAMP2"),label=T) #Endothelials
    FeaturePlot(dat,features=c("EPCAM","KRT19","CDH1","EGFR"),label=T) # Epithelials
    FeaturePlot(dat,features=c("TPSAB1","TPSB2","MS4A2","KIT"),label=T) # Mast cells
    FeaturePlot(dat,features=c("CLEC10A","CD1C","CLEC4C", "IL3RA"),label=T) # pDC 
    dev.off()

dat$infercnv.type <- "tmp"
dat$infercnv.type[which(dat$seurat_clusters%in%c(0,13,14))] <- "Tcell"
dat$infercnv.type[which(dat$seurat_clusters%in%c(7,9,11,17))] <- "Tumor"
saveRDS(dat,file="/public/workspace/lily/Mutiple_LB/E0927.RDS")




############################################################################################################################
# run infercnv 
library(infercnv)
library(Seurat)
E0927 <- readRDS("/public/workspace/lily/Mutiple_LB/E0927.RDS")
sub.dat <- subset(E0927,cells=which(E0927$infercnv.type%in%c("Tcell","Tumor")))


setwd('/public/workspace/lily/Mutiple_LB/infercnv/')
dir.create('./E0927')
write.table(sub.dat$infercnv.type,file="E0927_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub.dat@assays$RNA@counts)
write.table(count,file="E0927_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/public/workspace/lily/Mutiple_LB/infercnv/E0927_count_exp.txt",
         annotations_file="/public/workspace/lily/Mutiple_LB/infercnv/E0927_cell_info.txt",
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names="Tcell")


setwd('/public/workspace/lily/Mutiple_LB/infercnv/E0927')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=10, #big
                             no_plot=F ,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )




































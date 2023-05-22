
# 2021-12-18
# analysis about paired-LCBM
#===========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_BM.RDS")
# Paired_LCBM 
# get T cells and Epithelial to run infercnv 
sub.dat <- subset(dat,cells=which(dat$seurat_clusters%in%c(16,19,1,9,10)))
# saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_BM_infercnv.RDS")
setwd('/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/')
dir.create('./Pair_BM')
sub.dat$infercnv.type <- paste0("T",sub.dat$seurat_clusters)
sub.dat$infercnv.type[which(sub.dat$seurat_clusters%in%c(16,19))] <- "Tcell"
write.table(sub.dat$infercnv.type,file="Pair_BM_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub.dat@assays$RNA@counts)
write.table(count,file="Pair_BM_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_BM_count_exp.txt",
         annotations_file="/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_BM_cell_info.txt",
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names="Tcell")


setwd('/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_BM')
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





#==============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_Lung.RDS")
# Paired_Lung 
# get T cells and Epithelial to run infercnv 
sub.dat <- subset(dat,cells=which(dat$seurat_clusters%in%c(12,20,7,15,19,22)))
# saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_BM_infercnv.RDS")
setwd('/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/')
dir.create('./Pair_Lung')
sub.dat$infercnv.type <- paste0("T",sub.dat$seurat_clusters)
sub.dat$infercnv.type[which(sub.dat$seurat_clusters%in%c(12,20))] <- "Tcell"
write.table(sub.dat$infercnv.type,file="Pair_Lung_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub.dat@assays$RNA@counts)
write.table(count,file="Pair_Lung_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_Lung_count_exp.txt",
         annotations_file="/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_Lung_cell_info.txt",
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names="Tcell")


setwd('/public/workspace/lily/Lung2Brain/Pair_LCBM/inferCNV/Pair_Lung')
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












#=============================================================================================================================================
# DimPlot 
#=============================================================================================================================================
library(Seurat)
# LCBM data 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_BM.RDS")
dat$celltype <- "Undefined"
dat$celltype[which(dat$seurat_clusters%in%c(16,19))] <- "T&NK"
dat$celltype[which(dat$seurat_clusters%in%c(3))] <- "B cell"
dat$celltype[which(dat$seurat_clusters%in%c(4,7,12,15))] <- "Myeloid"
dat$celltype[which(dat$seurat_clusters%in%c(0,2,5,6,8))] <- "Fibro."
dat$celltype[which(dat$seurat_clusters%in%c(11))] <- "Endothelial"
dat$celltype[which(dat$seurat_clusters%in%c(9,10,1))] <- "Tumor"
dat$celltype[which(dat$seurat_clusters%in%c(18))] <- "Mast"


cols <- c('#377EB8','#910241','#984EA3',"#E7298A",'#F29403',"#B2DF8A",'#E41A1C','#999999')
DimPlot(dat,group.by="celltype",cols=cols)








# Lung data 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Pair_LCBM/Pair_Lung.RDS")
dat$celltype <- "Undefined"
dat$celltype[which(dat$seurat_clusters%in%c(12,20))] <- "T&NK"
dat$celltype[which(dat$seurat_clusters%in%c(14))] <- "B cell"
dat$celltype[which(dat$seurat_clusters%in%c(0,1,2,3,4,5,6,7,8,9,10,11,18))] <- "Myeloid"
dat$celltype[which(dat$seurat_clusters%in%c(13,16,19))] <- "Fibro."
dat$celltype[which(dat$seurat_clusters%in%c(21))] <- "Endothelial"
dat$celltype[which(dat$seurat_clusters%in%c(7,15))] <- "Tumor"
dat$celltype[which(dat$seurat_clusters%in%c(19,22))] <- "Epithelial"
dat$celltype[which(dat$seurat_clusters%in%c(17))] <- "Mast"


cols <- c('#377EB8','#910241',"#FB9A99",'#984EA3',"#E7298A",'#F29403',"#B2DF8A",'#E41A1C')
DimPlot(dat,group.by="celltype",cols=cols)




























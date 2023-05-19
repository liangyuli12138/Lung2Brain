
# 2021-12-27
# this program is used to run inferCNV for each samples.
# and then to identify Tumor cells
#============================================================================================================================================
# bytlib load JAGS-4.3.0
# R

library(Seurat)
library(infercnv)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")
dat$celltype <- "unknow"
dat$celltype[which(dat$seurat_clusters%in%c(0,1,21,23))] <- "Tcell"
dat$celltype[which(dat$seurat_clusters%in%c(3,6,8,11,12,14,19,22))] <- "Myeloid"
dat$celltype[which(dat$seurat_clusters%in%c(15,18))] <- "Bcell"
dat$celltype[which(dat$seurat_clusters%in%c(7,10,16,9))] <- "Fibroblast"
dat$celltype[which(dat$seurat_clusters%in%c(17))] <- "Oligo."
dat$celltype[which(dat$seurat_clusters%in%c(20))] <- "Endothelial"
dat$celltype[which(dat$seurat_clusters%in%c(2,5,13,4))] <- "Epithelial"


subdat <- subset(dat,cells=which(dat$celltype %in% c("Tcell","Epithelial")))


samplelist <- unique(subdat$orig.ident)
for(i in 1:length(samplelist)){

    tmp.dat <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
    DefaultAssay(tmp.dat) <- "RNA"
    dir.create(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i]))
    setwd(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i]))
    # respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])
    samplename <- samplelist[i]

    write.table(tmp.dat$celltype,paste(samplename,"_cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
    count<-as.matrix(tmp.dat@assays$RNA@counts)
    write.table(count,paste(samplename,"_count_exp.txt",sep='_'),sep="\t",quote=F)

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(samplename,"_count_exp.txt",sep='_'),
            annotations_file=paste(samplename,"_cell_info.txt",sep='_'),
            delim="\t",
            gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
            ref_group_names=c("Tcell"))
            
    infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="./", 
                                cluster_by_groups=T, 
                                plot_steps=F,
                                no_prelim_plot = TRUE,
                                num_threads=8, #big
                                no_plot=F ,
                                output_format = "pdf" # maybe can more quick 
                                # used for final scaling to fit range (0,2) centered at 1.
                                )

}





#============================================================================================================================================
# Epithelial re-inte and re-cluster
# 2021-12-27
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_V6.RDS")
subdat <- subset(dat,cells=which(dat$celltype=="Epithelial"))

inte.list <- list()
samplelist <- unique(subdat$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter = 50)
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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")





#===========================================================================================================================================
# 2021-12-27 
# to calculate CNV information and identify tumor cells 
#===========================================================================================================================================
library(Seurat)
library(infercnv)

samplelist<- c("A20190305","A20190312","BT1291","BT1292","BT1296","BT1297","D0927","E0927","Pair_BM","Pair_LUNG","scrBT1431m","scrBT1432m","T_Bsc1")
tmp.list <- list()
for(i in 1:length(samplelist)){
    respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i],"/")
    dat <- readRDS(paste0(respath,"run.final.infercnv_obj"))
    # calculate all cells CNV.score
    tmp.cnv <- apply(dat@expr.data,2,function(x){sum((x-1)^2)})

    # get reference CNV.score
    # cutoff use 75% of reference CNV score 
    ref.cnv <- quantile(tmp.cnv[dat@reference_grouped_cell_indices$Tcell],0.75)

    # get obs cells 
    obs <- tmp.cnv[dat@observation_grouped_cell_indices$Epithelial]
    obs.res <- obs[which(obs>ref.cnv)]
    # print(samplelist[i])
    # print(length(dat@observation_grouped_cell_indices$Epithelial))
    # print(length(obs.res))
    # tmp.list[[i]] <- obs.res
    # names(tmp.list)[i] <- samplelist[i]
    saveRDS(obs.res,file=paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i],"_tumorcell_cnv.RDS"))
}






#==============================================================================================================================================
# 2021-12-27 
# use CNV information to identify Tumor cells 
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")

samplelist<- c("A20190305","A20190312","BT1291","BT1292","BT1296","BT1297","D0927","E0927","Pair_BM","Pair_LUNG","scrBT1431m","scrBT1432m","T_Bsc1")
tmp.list <- list()
for(i in 1:length(samplelist)){
    tmp.cnv <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i],"_tumorcell_cnv.RDS"))
    tmp.list[[i]] <- tmp.cnv
    names(tmp.list)[i] <- samplelist[i]
}

# get cell name 
cell.tumor <- unlist(tmp.list)
names(cell.tumor) <- sapply(strsplit(names(cell.tumor),"\\."),function(x){x[[2]]})

# annoatation
dat$celltype.refine <- "Epithelial_like"
dat$celltype.refine[which(colnames(dat)%in%names(cell.tumor))] <- "Tumor"

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")




















#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
# 2021-12-30 

library(Seurat)
library(infercnv)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
subdat <- subset(dat,cells=which(dat$celltype %in% c("Tcell","Epithelial")))

samplelist <- unique(subdat$orig.ident)
for(i in 1:length(samplelist)){

    tmp.dat <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
    DefaultAssay(tmp.dat) <- "RNA"
    dir.create(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i]))
    setwd(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i]))
    # respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])
    samplename <- samplelist[i]

    write.table(tmp.dat$celltype,paste(samplename,"_cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
    count<-as.matrix(tmp.dat@assays$RNA@counts)
    write.table(count,paste(samplename,"_count_exp.txt",sep='_'),sep="\t",quote=F)

    infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(samplename,"_count_exp.txt",sep='_'),
            annotations_file=paste(samplename,"_cell_info.txt",sep='_'),
            delim="\t",
            gene_order_file="/public/workspace/lily/REF/INDEX-hg38/hg38_position_pure.txt",
            ref_group_names=c("Tcell"))
            
    infercnv_obj = infercnv::run(infercnv_obj,
                                cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                out_dir="./", 
                                cluster_by_groups=T, 
                                plot_steps=F,
                                no_prelim_plot = TRUE,
                                num_threads=8, #big
                                no_plot=F ,
                                output_format = "pdf" # maybe can more quick 
                                # used for final scaling to fit range (0,2) centered at 1.
                                )

}




#===========================================================================================================================================
# 2022-1-1 
# re-inte and re cluster Epithelials 
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
subdat <- subset(tmp.dat,cells=which(tmp.dat$celltype=="Epithelial"))

inte.list <- list()
samplelist <- unique(subdat$orig.ident)
for(i in 1: length(samplelist)){
    tmp <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter = 30)
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

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")





#=============================================================================================================================================
# 2022-1-1 
# move hg19 result into a hg19 files
# use CNV information to identify Tumor cells 
library(Seurat)
library(infercnv)
samplelist<- gsub("\\.RDS$","",grep("*.RDS$",dir("/public/workspace/lily/Lung2Brain/Version6/Prepare_Data/"),value=T))
tmp.list <- list()
for(i in 1:length(samplelist)){
    respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i],"/")
    dat <- readRDS(paste0(respath,"run.final.infercnv_obj"))
    # calculate all cells CNV.score
    tmp.cnv <- apply(dat@expr.data,2,function(x){sum((x-1)^2)})

    # get reference CNV.score
    # cutoff use 75% of reference CNV score 
    ref.cnv <- quantile(tmp.cnv[dat@reference_grouped_cell_indices$Tcell],0.75)

    # get obs cells 
    obs <- tmp.cnv[dat@observation_grouped_cell_indices$Epithelial]
    obs.res <- obs[which(obs>ref.cnv)]
    print(samplelist[i])
    print(length(dat@observation_grouped_cell_indices$Epithelial))
    print(length(obs.res))
    tmp.list[[i]] <- obs.res
    names(tmp.list)[i] <- samplelist[i]
    saveRDS(obs.res,file=paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i],"_tumorcell_cnv.RDS"))
}


#=============================================================================================================================================
# and now to define cells 
#=============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
samplelist<- gsub("\\.RDS$","",grep("*.RDS$",dir("/public/workspace/lily/Lung2Brain/Version6/Prepare_Data/"),value=T))
tmp.list <- list()
for(i in 1:length(samplelist)){
    tmp.cnv <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i],"_tumorcell_cnv.RDS"))
    tmp.list[[i]] <- tmp.cnv
    names(tmp.list)[i] <- samplelist[i]
}

# get cell name 
cell.tumor <- unlist(tmp.list)
names(cell.tumor) <- sapply(strsplit(names(cell.tumor),"\\."),function(x){x[[2]]})

# annoatation
dat$celltype.refine <- "Epithelial_like"
dat$celltype.refine[which(colnames(dat)%in%names(cell.tumor))] <- "Tumor"

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")






# 2022-2-22
# change inte16s sample to add tumor and non-tumor information
#===================================================================================================================================

library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
epi <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")

tmp.dat$celltype.refine <- tmp.dat$celltype
tmp.dat$celltype.refine[which(colnames(tmp.dat)%in%names(which(epi$celltype.refine=="Tumor")))] <- "Tumor"

saveRDS(tmp.dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")







































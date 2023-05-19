# 2022-6-15
# Figure 0

# this program is used to analysis paired LCBM data
# 2022-5-6 try to use scenic result to find some consequencent result 
# 1. velocity result in Analysis_velocyto.R 
# 2. integration paired sample (for tumor cells , for all non-Epi just get subset)
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
inte.list<- list()
samplelist <- c("Pair_BM","Pair_Lung")
for(i in 1:length(samplelist)){
  tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]&dat$celltype.refine=="Tumor"))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
# now integration 
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
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")


# and plot a cell type plot 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_V6","HPSC_C5"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

dat$BMS_score <- mod$BMS_V6_norm
tmp <- aggregate(BMS_score~seurat_clusters,data=dat@meta.data,FUN=median)
quantile(tmp$BMS_score) # use 75% and 25% as cutoff 4,5,6 is BMS high and 0,8,9 is BMS low 
dat$BMS.type <- "MID"
dat$BMS.type[which(dat$seurat_clusters%in%c(4,5,6))] <- "High"
dat$BMS.type[which(dat$seurat_clusters%in%c(0,8,9))] <- "Low"

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")


# plot result 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/Paired_tumor_BMSscore_type_group.pdf",useDingbats=F)
VlnPlot(dat,features="BMS_score",group.by="type_group",pt.size=0)

ggplot(data=dat@meta.data,aes(x=type_group,y=BMS_score,fill=type_group))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.3, color="black",outlier.colour = NA) 

dev.off()


# analysis non-tumor cell type change 
# so many stroma
#==========================================================================================================================================
#==========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
subdat <- subset(dat,cells=which(dat$orig.ident%in%c("Pair_BM","Pair_Lung")))

##### insert start
# 2022-6-27
# re-integration these sample
inte.list <- list()
samplelist <- unique(subdat$orig.ident)
for(i in 1:length(samplelist)){
    tmp <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
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
inte <- FindClusters(inte,resolution=1)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version6/Data/Pair_Sample_inte.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_DimPlot.pdf",useDingbats=F,width=10,height=8)
cols <- c('#377EB8','#910241','#fc636b','#984EA3','#F29403','#aea400','#B2DF8A')
DimPlot(inte,group.by="celltype",split.by="type_group",cols=cols,raster=F)
dev.off()
#### insert done



dat.res <- subset(subdat,cells=colnames(subdat)[-which(subdat$celltype=="Epithelial")]) # do not use Epithelial
# analysis percentage differernce
dat <- subset(dat.res,cells=which(dat.res$celltype%in%c("Bcell","Myeloid","Oligo.","Tcell")))
tmp <- table(dat$type_group,dat$celltype)
res.f <- apply(tmp,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Myeloid","Oligo.","Tcell"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/Immunecells_celltype.pdf",useDingbats=F)
cols <- c('#377EB8','#F29403','#aea400','#B2DF8A')
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()




#=============================================================================================================================================
dat.res <- subset(subdat,cells=colnames(subdat)[-which(subdat$celltype=="Epithelial")]) # do not use Epithelial
# analysis percentage differernce
dat <- dat.res
tmp <- table(dat$type_group,dat$celltype)
res.f <- apply(tmp,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Endothelial","Fibroblast","Myeloid","Oligo.","Tcell"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/non_Epithelial_celltype.pdf",useDingbats=F)
cols <- c('#377EB8','#910241','#984EA3','#F29403','#aea400','#B2DF8A')
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()








# analysis scenic result ,try to find some consequencent result 
# 2022-5-6
# 2022-6-10 lly decided not use PySCENIC result ,use dorothea 
#==========================================================================================================================================
#==========================================================================================================================================
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
# tmp.dat <- read.csv("/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/step3.auc_mtx.csv")

# rownames(tmp.dat) <- tmp.dat$Cell
# tmp.dat$Cell <- NULL
# colnames(tmp.dat) <- gsub("\\.\\.\\.$","",colnames(tmp.dat))
# # check cell names
# all(colnames(dat)==rownames(tmp.dat))
# tmp.dat$type_group <- dat$type_group

# # calculate by group
# tmp.res <- data.frame(aggregate(.~type_group,data=tmp.dat,FUN=median))
# rownames(tmp.res) <- tmp.res$type_group
# tmp.res$type_group <- NULL
# tmp.res <- data.frame(t(tmp.res))
# final.res <- data.frame(apply(tmp.res,2,function(x){as.numeric(as.vector(x))}))
# rownames(final.res) <- rownames(tmp.res)



# # do some filter

# res.filter <- final.res[which(final.res$LCBM > final.res$nMLUAD & final.res$MLUAD > final.res$nMLUAD),]
# res.filter$LCBM.nmLUAD <- res.filter$LCBM - res.filter$nMLUAD
# res.filter$MLUAD.nmLUAD <- res.filter$MLUAD - res.filter$nMLUAD









# 2022-5-9 
# use another ways to calculate TFs
# need bytlib load r/4.0.4
# 2022-6-8
#==========================================================================================================================================
#==========================================================================================================================================

library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

# dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")

# dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))
# regulon <- dorothea_regulon_human %>% dplyr::filter(confidence %in% c("A","B","C"))
# # check DefaultyAssay 
# dat <- run_viper(dat, regulon, assay_key = "RNA",options = list(method = "scale", minsize = 4, 
#                                  eset.filter = FALSE, cores = 20, 
#                                  verbose = FALSE))


# res <- data.frame(t(dat@assays$dorothea@data))
# res$group <- dat$type_group
# saveRDS(res,file="/public/workspace/lily/Lung2Brain/Version6/Data/Fig3/Dorothe_res_Tumor_scRNA.RDS")
res <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Fig3/Dorothe_res_Tumor_scRNA.RDS")

aggregate(.~group,data=res,FUN=median) -> tmp.res
rownames(tmp.res) <- tmp.res[,1]
tmp.res <- tmp.res[,-1]
tmp.res <- data.frame(t(tmp.res))

res.filter <- tmp.res[which(tmp.res$LCBM > tmp.res$nMLUAD & tmp.res$MLUAD > tmp.res$nMLUAD),]
res.filter$LCBM.nmLUAD <- res.filter$LCBM - res.filter$nMLUAD
res.filter$MLUAD.nmLUAD <- res.filter$MLUAD - res.filter$nMLUAD



# now use rabitt data to do filter
rabit.dat <- read.table("/public/workspace/lily/metastasis/data/TCGA_RABIT/LUAD/RABIT_LUAD.HiSeq.V2",header=T,sep="\t")
rownames(rabit.dat) <- gsub("-",".",rabit.dat[,1])
rabit.dat$X <- NULL
rabit.dat <- as.matrix(rabit.dat)
# load BMS score
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
rabit.dat.f <- rabit.dat[intersect(rownames(rabit.dat),rownames(mod)),]
mod.f <- mod[intersect(rownames(rabit.dat),rownames(mod)),]

tmp.res <- t(apply(rabit.dat.f,2,function(x){
  c(cor.test(as.numeric(x),mod.f$BMS_V6_norm,method="spearman")$estimate,
  cor.test(as.numeric(x),mod.f$BMS_V6_norm,method="spearman")$p.value)
  }))

tmp.res.f <- tmp.res[which(tmp.res[,2]<0.05),]
colnames(tmp.res.f) <- c("Rho","Pvalue")

# now combine dorothea data 
res.dat <- merge(res.filter,tmp.res.f,by="row.names")
res.dat$group <- "negative"
res.dat$group[which(res.dat$Rho>0)] <- "positive"


# plot result 
library(ggplot2)
library(ggrepel)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/TFs_Tumor.pdf",useDingbats=F)
ggplot(res.dat,aes(x=LCBM.nmLUAD,y=Rho,color=group)) + geom_point(aes(size=MLUAD.nmLUAD)) + 
    scale_color_manual(values=c("#2e9df7","#ec2c22")) + 
    geom_text(aes(label = Row.names), size = 3)+ labs(title="TFs")
dev.off()

















# 2022-6-9
# add scvelo result
# calculate in Paired data
# 2022-5-10 monocle show not good result 
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res_time.RDS")
dat$cellname <- sapply(strsplit(colnames(dat),"_"),function(x){x[[1]]})
rownames(tmp.dat)<- tmp.dat$cellname
tmp.dat <- tmp.dat[dat$cellname,]
all(dat$cellname==rownames(tmp.dat))
dat$pseudotime <- tmp.dat$pseudotime
dat$latent_time <- tmp.dat$latent_time
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")


# 1. plot a barplot for this pseudotime
library(ggplot2)
plot.dat <- dat@meta.data[,c("cellname","type_group","pseudotime")]
plot.dat <- plot.dat[order(plot.dat$pseudotime),]
plot.dat$cellname <- factor(plot.dat$cellname,levels=plot.dat$cellname)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/Pair_Tumor_scvelo.pdf",useDingbats=F)
ggplot(data=plot.dat,aes(x=cellname,y=pseudotime,color=type_group))+geom_bar(stat="identity",width=0)
dev.off()



# 1.2 calculate cytoTrace
library(Seurat)
library(CytoTRACE)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")

# cell stemness cytoTrace
results <- CytoTRACE(as.matrix(dat[["RNA"]]@data),ncores = 4,subsamplesize = 1000)
dat$cytotrace <- results$CytoTRACE

plot.dat <- dat@meta.data[,c("cellname","type_group","pseudotime","cytotrace")]
plot.dat <- plot.dat[order(plot.dat$pseudotime),]
plot.dat$cellname <- factor(plot.dat$cellname,levels=plot.dat$cellname)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/Pair_Tumor_scvelo.pdf",useDingbats=F)
ggplot(data=plot.dat,aes(x=cellname,y=pseudotime,color=type_group))+geom_bar(stat="identity",width=0)
dev.off()




# just plot density 
# not OK
# ggplot(plot.dat,aes(x=pseudotime))+geom_density(aes(color=type_group))






# 2. use GSE200563 data to calculate TIP signature 
# not good 
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("MLUNG","BM")),]
dat <- tmp.dat[,rownames(sampleinfo)]

# calculate 
modname <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/TIP/"))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),modname,"/public/workspace/lily/MOD_file/TIP/",permN=0)
mod <- as.data.frame(mod)




source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
tmp.res <- readRDS("~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")
modname <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/TIP/"))
mod.gse <- as.data.frame(mod.analyze2(as.matrix(tmp.res$GSE_Data),modname,"/public/workspace/lily/MOD_file/TIP/",permN=0))
mod.tcga <- as.data.frame(mod.analyze2(as.matrix(tmp.res$TCGA_Data),modname,"/public/workspace/lily/MOD_file/TIP/",permN=0))


mod.gse$group <- tmp.res$GSE_info$group
tmp.1 <- aggregate(.~group,data=mod.gse[,c(26:50,76)],FUN=median)
tmp.2 <- apply(mod.tcga[,26:50],2,function(x){median(x)})
all(names(tmp.1)[-1]==names(tmp.2))

# merge result 
rownames(tmp.1) <- tmp.1[,1]
tmp.1 <- tmp.1[,-1]
tmp.res <- data.frame(t(tmp.1))
tmp.res$TCGA <- unname(tmp.2)









# 3. calculate Metabolism pathway in TCGA MLUG and BM
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
tmp.res <- readRDS("~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")
modname <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/"))
mod.gse <- as.data.frame(mod.analyze2(as.matrix(tmp.res$GSE_Data),modname,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0))
mod.tcga <- as.data.frame(mod.analyze2(as.matrix(tmp.res$TCGA_Data),modname,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0))

mod.gse$group <- tmp.res$GSE_info$group
tmp.1 <- aggregate(.~group,data=mod.gse[,c(86:170,256)],FUN=median)
tmp.2 <- apply(mod.tcga[,86:170],2,function(x){median(x)})
all(names(tmp.1)[-1]==names(tmp.2))

# merge result 
rownames(tmp.1) <- tmp.1[,1]
tmp.1 <- tmp.1[,-1]
tmp.res <- data.frame(t(tmp.1))
tmp.res$TCGA <- unname(tmp.2)






















#============================================================================================================================
# 2022-6-25
# plot trajectory 

library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_Tumor_1.5K_monocle.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")

all(colnames(dat)==colnames(tmp.dat))
dat$Pseudotime <- tmp.dat$pseudotime

library(future)
future::plan(multisession, workers=10)
diff_test_res <- differentialGeneTest(dat,fullModelFormulaStr = "~sm.ns(Pseudotime)")
tmp.res <- diff_test_res[which(diff_test_res$qval<0.05),]
tmp.res <- tmp.res[order(tmp.res$qval,decreasing=F),]

p1 <- plot_pseudotime_heatmap(dat[rownames(tmp.res)[1:50],],cluster_rows=T,num_clusters =3, cores = 5,show_rownames = T,return_heatmap=T)







# 2022-6-27
# calculate hallmark
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")

modname <- gsub("\\.mod$","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),modname,"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)

tmp.res <- mod[,51:100]
tmp.res$group <- dat$type_group
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Paired_Tumor_hallmark_mod.RDS")

HALLMARK_ANGIOGENESIS_norm
HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm

library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Hallmark_H_LCBM_Pair.pdf",useDingbats=F)
ggplot(tmp.res,aes(x=HALLMARK_ANGIOGENESIS_norm))+geom_density(aes(fill=group),alpha=0.4)
ggplot(tmp.res,aes(x=HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm))+geom_density(aes(fill=group),alpha=0.4)
dev.off()




# metabolism
modname <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/"))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),modname,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod <- as.data.frame(mod)

tmp.res <- mod[,86:170]
tmp.res$group <- dat$type_group
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Paired_Tumor_metabolism_mod.RDS")

pdf("~/tmp.pdf")
apply(tmp.res[,-15],2,function(x){
  p <- sm.density.compare(as.numeric(x),factor(tmp.res$group))
  print(p)
})
dev.off()














#========================================================================================================================
# plot genome stability in Paired Samples
library(infercnv)
library(Seurat)

tmp.dat <- readRDS("/public/workspace/lily/metastasis/data/hg38_gene_Chrpq.RDS")
tumor <- readRDS("~/Lung2Brain/Version6/Data/Tumor.RDS")
# change function
scCNV2chr <- function(sampleid){
    dat <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",sampleid,"/run.final.infercnv_obj"))
    expr <- as.data.frame(dat@expr.data[,-dat@reference_grouped_cell_indices$Tcell])
    expr <- expr[,which(colnames(expr)%in%colnames(tumor))]
    tmp.dat <- tmp.dat[rownames(expr),]
    expr$gene_order <- tmp.dat$chrpq
    res <- aggregate(.~gene_order,data=expr,FUN=mean) #get every chr for each cell 
    rownames(res) <- res$gene_order
    res$gene_order <- NULL
    # res.f <- apply(res,1,function(x){mean(x)})
    return(res)
}

scCNV.PLung <- scCNV2chr("Pair_Lung")

scCNV.PBM <- scCNV2chr("Pair_BM")

tmp.res <- data.frame(
  MADCNV=as.numeric(c(apply(scCNV.PLung,2,function(x){mad(x)}),apply(scCNV.PBM,2,function(x){mad(x)}))),
  group=factor(c(rep("PLung",ncol(scCNV.PLung)),rep("PBM",ncol(scCNV.PBM))))
  )

library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_MADCNV.pdf",useDingbats=F)
ggplot(data=tmp.res,aes(x=group,y=MADCNV,fill=group))+
  geom_violin(width=1.4) +
  geom_boxplot(width=0.3, color="black",outlier.colour = NA) 

dev.off()









#========================================================================================================================
# plot TCGA sample 
# use TCGA early samples
#===========================================================================================================================
# 2022-6-22
# analysis all chr gene signature in paired tumor cells
# bytlib load r/4.1.2
# library(msigdb)
# geneset <- msigdb::getMsigdb(org="hs",id="SYM")
# # names(geneset)[1:261]
# # all chr gene 
# genelist <- matrix(ncol=2)
# for(i in 1:261){
#   genelist <- rbind(genelist,c(geneset[[i]]@setName,paste(geneset[[i]]@geneIds,collapse=";")))
# }
# # get gene signature
# genelist <- data.frame(genelist)
# genelist <- genelist[-1,]
# colnames(genelist) <- c("setName","gene")
# # integration
# genelist$group <- stringr::str_extract(genelist$setName,"chr[0-9]+[pq]")
# tmp.res <- (aggregate(gene~group,data=genelist,function(x){paste(x,collapse=";")}))
# saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/MSGDB_chr_gene.RDS")

# create mod
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')

for(i in 1:nrow(tmp.res)){
  gene <- unique(unlist(strsplit(tmp.res[i,2],";")))
  mod.generate(gene,tmp.res[i,1],out=paste0("/public/workspace/lily/MOD_file/MSGDB_Chr/",tmp.res[i,1],".mod")) # make a mod file 
}


#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
#===========================================================================================================================
# calculate chr gene 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_sample_tumor_inte.RDS")
modname <- gsub("\\.mod$","",dir("/public/workspace/lily/MOD_file/MSGDB_Chr/"))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),modname,"/public/workspace/lily/MOD_file/MSGDB_Chr/",permN=0)
mod <- as.data.frame(mod)
tmp.res <- mod[,41:80]
PTCNV <- tail(sort(apply(tmp.res,2,function(x){mean(x)})),10)
CNVdata <- as.data.frame(PTCNV)
CNVdata$chr <- sapply(strsplit(rownames(CNVdata),"_"),function(x){x[[1]]})
CNVdata$chr <- factor(CNVdata$chr,levels=c(CNVdata$chr))

label_data <- CNVdata
# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$PTCNV-0.5) /number_of_bar  
# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_Chr_Top10.pdf",useDingbats=F)
ggplot(CNVdata, aes(x=chr, y=PTCNV)) +    
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  
  #The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-0.2,1.5) +
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")     # This remove unnecessary margin around plot
  ) +
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0)+
  geom_text(data=label_data, aes(x=chr, y=PTCNV+0.1, label=chr, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

dev.off()




# use Top10 and CNV info 
# 增加了CNV的信息，用这个排序
# 既考虑了在转移的肿瘤细胞中都要高，又考虑了要在脑转移灶更高
tmpdat <- data.frame(PBM=as.numeric(apply(scCNV.PBM,1,function(x){median(x)})),PLung=as.numeric(apply(scCNV.PLung,1,function(x){median(x)})))
rownames(tmpdat) <- names(apply(scCNV.PBM,1,function(x){median(x)}))
tmpdat$diff <- tmpdat$PBM - tmpdat$PLung
tmpdat.f <- tmpdat[which(rownames(tmpdat)%in%c("chr1p","chr1q","chr2p","chr2q","chr3p","chr3q","chr4p","chr4q","chr5p","chr5q")),]
tmpdat.f <- tmpdat.f[order(tmpdat.f$diff,decreasing=T),]
tmpdat.f$chr <- factor(rownames(tmpdat.f),levels=rownames(tmpdat.f))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_Chr_Top10_CNV.pdf",useDingbats=F)
ggplot(tmpdat.f, aes(x=chr, y=diff)) +    
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("blue", 0.3)) +
  coord_flip() + theme_bw()

dev.off()

















############################################################################################################################
# calculate in TCGA early stage samples
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
tmp.res <- readRDS("~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")
modname <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/MSGDB_Chr/"))
mod.gse <- as.data.frame(mod.analyze2(as.matrix(tmp.res$GSE_Data),modname,"/public/workspace/lily/MOD_file/MSGDB_Chr/",permN=0))
mod.tcga <- as.data.frame(mod.analyze2(as.matrix(tmp.res$TCGA_Data),modname,"/public/workspace/lily/MOD_file/MSGDB_Chr/",permN=0))

mod.gse$group <- tmp.res$GSE_info$group
tmp.tcga <- mod.tcga[,41:80]
tmp.gse <- mod.gse[,c(41:80,121)]

tmp.res <- aggregate(.~group,data=tmp.gse,FUN=median)
rownames(tmp.res) <- tmp.res$group
tmp.res$group <- NULL
tmp.res <- data.frame(t(tmp.res))
all(rownames(tmp.res)==names(apply(tmp.tcga,2,function(x){median(x)})))
# get TCGA data 
tmp.res$TCGA <- apply(tmp.tcga,2,function(x){median(x)})
saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/TCGA_GSE200563_Chr_mod.RDS")

#===============================================================================================================
# plot result
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/TCGA_GSE200563_Chr_mod.RDS") 
tmp.res <- tmp.res[order(tmp.res$TCGA),]
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/TCGA_early_GSE200563_Chr.pdf",useDingbats=F)
plot(tmp.res$MLUNG,tmp.res$TCGA)
abline(h=0.4)
abline(v=0.4)
dev.off()


















#===========================================================================================================================
# calculate T cells activation signature 
# 2022-6-28
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_Sample_inte.RDS")
Tsub <- subset(dat,cells=which(dat$celltype.refine=="Tcell"))
RE_INTE <- function(dat,sample){
    inte.list <- list()
    samplelist <- unique(dat@meta.data[,sample])
    for(i in 1:length(samplelist)){
        tmp <- subset(dat,cells=which(dat@meta.data[,sample]==samplelist[i]))
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
    inte <- FindClusters(inte,resolution=1)
    #TSNE
    # if Umap can not use
    inte <- RunTSNE(inte)

    return(inte)
}
Tdat <- RE_INTE(Tsub,sample="orig.ident")
Tdat$celltype.refine[which(Tdat$seurat_clusters%in%c(2,3,5))] <- "NK"
saveRDS(Tdat,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Pair_Sample_Tdat.RDS")
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_T_NK_DimPlot.pdf",useDingbats=F)
DimPlot(Tdat,group.by="celltype.refine")
dev.off()


Tsub <- subset(Tdat,cells=which(Tdat$seurat_clusters%in%c(0,1,4,6,7,8)))

Tact <- c("STAT1","GNLY","GZMA","GZMB","GZMK","GZMM","CD70","CD27",
"IRF1","IRF4","IRF7","IFI27","IFI35","NFATC2","NFATC1","CD28","CD40LG",
"MX1","STAT3","NT5E","HAVCR2","TNFRSF18","LAG3","IRF2","IRF3","IRF5","IRF8",
"JAK2","NFATC4","NFATC3","CD24","CD44")
Tdys <- c("PDCD1","HAVCR2","TGFB1","TGFB2","TGFB3","TIGIT","CTLA4")


# make gene mod 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(Tact,"T_act",out="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/Tcell_act.mod") # make a mod file 
mod.generate(Tdys,"T_dys",out="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/Tcell_dys.mod") # make a mod file 

# calculate score 

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(Tsub[['RNA']]@data),c("Tcell_act","Tcell_dys"),"/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$group <- Tsub$type_group
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Pair_Sample_Tsub_act_dys_mod.RDS")


# plot result 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_Tsub_act_dys_boxplot.pdf",useDingbats=F)
boxplot(Tcell_act_norm~group,data=mod,FUN=median,outline=F,ylim=c(0,1))
boxplot(Tcell_dys_norm~group,data=mod,FUN=median,outline=F,ylim=c(0,1))
dev.off()











#===========================================================================================================================
# calculate Myeloid cells activation signature 
# 2022-6-28
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Pair_Sample_inte.RDS")
Msub <- subset(dat,cells=which(dat$celltype.refine=="Myeloid"))
RE_INTE <- function(dat,sample){
    inte.list <- list()
    samplelist <- unique(dat@meta.data[,sample])
    for(i in 1:length(samplelist)){
        tmp <- subset(dat,cells=which(dat@meta.data[,sample]==samplelist[i]))
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
    inte <- FindClusters(inte,resolution=1)
    #TSNE
    # if Umap can not use
    inte <- RunTSNE(inte)

    return(inte)
}
Mdat <- RE_INTE(Msub,sample="orig.ident")
DefaultAssay(Mdat) <- "RNA"
Mdat$celltype.refine <- "Macrophage"
Mdat$celltype.refine[which(Mdat$seurat_clusters%in%c(15))] <- "Mast"
# FeaturePlot(Mdat,features=c("CD68","KIT","PPBP","LYZ","CST3","CD163","CSF3R"))
saveRDS(Mdat,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Pair_Sample_Mdat.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_Myeloid_DimPlot.pdf",useDingbats=F)
DimPlot(Mdat,group.by="celltype.refine")
dev.off()


Msub <- subset(Mdat,cells=which(Mdat$seurat_clusters%in%c(0:14,16:18)))

M1 <- c("IL6","IL12A","TNF","CD80","CD86")
M2 <- c("ARG1","MRC1","CCL17","CCL22")


# make gene mod 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(M1,"M1",out="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/M1.mod") # make a mod file 
mod.generate(M2,"M2",out="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/M2.mod") # make a mod file 

# calculate score 

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(Msub[['RNA']]@data),c("M1","M2"),"/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$group <- Msub$type_group
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Pair_Sample_Msub_M1M2_mod.RDS")


# plot result 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/Pair_Sample_Msub_M1M2_boxplot.pdf",useDingbats=F)
boxplot(M1_norm~group,data=mod,FUN=median,outline=F,ylim=c(0,1))
boxplot(M2_norm~group,data=mod,FUN=median,outline=F,ylim=c(0,1))
dev.off()







#===========================================================================================================================================
# GSE200563 data verify 
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("TIME_BM","TIME_LUNG","TME_BM") ),]
# & tmp.sampleinfo$Patient %in%c("Patient_15","Patient_19","Patient_35","Patient_5")
dat <- tmp.dat[,rownames(sampleinfo)]

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("M1","M2","Tcell_act","Tcell_dys"),"/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$group <- sampleinfo$group
mod$group <- as.vector(mod$group)
mod$group[grep("BM",mod$group)] <- "BM"

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/GSE200563_sig_verify.pdf",useDingbats=F)
boxplot(M1_norm~group,data=mod,FUN=median)
boxplot(M2_norm~group,data=mod,FUN=median)
boxplot(Tcell_dys_norm~group,data=mod,FUN=median)
boxplot(Tcell_act_norm~group,data=mod,FUN=median)
dev.off()


# mins
mod$diffT <- mod$Tcell_act_norm -mod$Tcell_dys_norm
# pvalue=0.11
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig0/GSE200563_Tminussig_verify.pdf",useDingbats=F)
boxplot(diffT~group,data=mod,FUN=median)
beeswarm::beeswarm(diffT~group,data=mod,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

dev.off()












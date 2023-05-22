
# this program is used to verify three subtype tumor cell is exist 
# 2021-7-2 
# 2021-11-23 use BMS to get high BMS and low BMS tumor cells
# 1. use a self data D and E 
# 2. use GSE 123904 data to verify 
#=============================================================================================================

# 1. use a self data to verify 
# D0927
library(Seurat) 
# dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
# calculate signature 
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)
# dat$BMS.update <- mod[,4]
# dat$Brain.gene <- mod[,5]
# dat$Lung.gene <- mod[,6]
# saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")

# library(ggplot2)
# pdf("/public/workspace/lily/Lung2Brain/Version5/Fig3/Fig3_D0927_Tumor.pdf",useDingbats=F)
# DimPlot(dat)
# FeaturePlot(dat,features="BMS.update",order=T)+
#     scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.7,0.7,1.0))  
# FeaturePlot(dat,features="Lung.gene",order=T) + 
#     scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.6,0.6,1.0))  
# FeaturePlot(dat,features="Brain.gene",order=T) +
#     scale_colour_gradientn(colours=c("steelblue","steelblue","white","#ff3c41","red"),values=c(0,0.45,0.45,1.0)) 
# # DimPlot(dat,reduction="tsne")
# dev.off()


####################################################################################################
# 2021-7-3 recluster then plot result 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")

recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	tmp_dat <- FindClusters(object = tmp_dat,resolution=1)
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	return(tmp_dat)
}
dat <- recluster(dat)

# define cell sub group 
tmp.group <- apply(dat@meta.data[,c("Brain.gene","Lung.gene","BMS.update")],1,function(x){
    tmp <- c()
    res <- c()
    if(x[1]>quantile(dat$Brain.gene,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"group16")
    }
    if(x[2]>quantile(dat$Lung.gene,0.75)){
        tmp <- c(tmp,x[2])
        res <- c(res,"group09")
    }
    if(x[3]>quantile(dat$BMS.update,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"BMS.update")
    }

    if(length(res)>0){
        res.f <- res[which.max(tmp)]
    }else{
        res.f <- "Unclassify"
    }

    return(res.f)
})
dat$tmp.group <- tmp.group
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/D0927_tumor_subgroup.pdf",useDingbats=F)
DimPlot(dat,group.by="tmp.group",cols=c("#ffc719","#2baf2b","#00a5dc","#ced7df"),reduction="tsne")
dev.off()



# 2021-11-23
dat$BMS.group <- "inte"
dat$BMS.group[which(dat$BMS.update<quantile(dat$BMS.update,0.25))] <- "low"
dat$BMS.group[which(dat$BMS.update>quantile(dat$BMS.update,0.75))] <- "high"
DefaultAssay(dat) <- "RNA"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")















##############################################################################################################################
# another sample E0927
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.5)
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	return(tmp_dat)
}
dat <- recluster(dat)

# define cell sub group 
tmp.group <- apply(dat@meta.data[,c("Brain.gene","Lung.gene","BMS.update")],1,function(x){
    tmp <- c()
    res <- c()
    if(x[1]>quantile(dat$Brain.gene,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"group16")
    }
    if(x[2]>quantile(dat$Lung.gene,0.75)){
        tmp <- c(tmp,x[2])
        res <- c(res,"group09")
    }
    if(x[3]>quantile(dat$BMS.update,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"BMS.update")
    }

    if(length(res)>0){
        res.f <- res[which.max(tmp)]
    }else{
        res.f <- "Unclassify"
    }

    return(res.f)
})
dat$tmp.group <- tmp.group
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")





pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/E0927_tumor_subgroup.pdf",useDingbats=F)
DimPlot(dat,group.by="tmp.group",cols=c("#ffc719","#2baf2b","#00a5dc","#ced7df"),reduction="tsne")
dev.off()


# 2021-11-23
dat$BMS.group <- "inte"
dat$BMS.group[which(dat$BMS.update<quantile(dat$BMS.update,0.25))] <- "low"
dat$BMS.group[which(dat$BMS.update>quantile(dat$BMS.update,0.75))] <- "high"
DefaultAssay(dat) <- "RNA"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")













######################################################################################################################
# MSK data show result 
# use GSE123902.RDS which inteagrate primary normal and brain metastasis samples.
#=====================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")

sub.dat <- subset(dat,cells=which(dat$seurat_clusters==6&dat$group=="METASTASIS"))
inte.list <- list()
samples <- unique(sub.dat$sample)
for(i in 1:length(samples)){
    tmp <- subset(sub.dat,cells=which(sub.dat$sample==samples[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}

# re - integration 
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=50,k.score=50)
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
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")

# calculate BMS score and Brain/Lung signature
#===============================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
dat$BMS.update <- mod[,4]
dat$Brain.gene <- mod[,5]
dat$Lung.gene <- mod[,6]

# define cell sub group 
tmp.group <- apply(dat@meta.data[,c("Brain.gene","Lung.gene","BMS.update")],1,function(x){
    tmp <- c()
    res <- c()
    if(x[1]>quantile(dat$Brain.gene,0.8)){
        tmp <- c(tmp,x[1])
        res <- c(res,"group16")
    }
    if(x[2]>quantile(dat$Lung.gene,0.8)){
        tmp <- c(tmp,x[2])
        res <- c(res,"group09")
    }
    if(x[3]>quantile(dat$BMS.update,0.8)){
        tmp <- c(tmp,x[1])
        res <- c(res,"BMS.update")
    }

    if(length(res)>0){
        res.f <- res[which.max(tmp)]
    }else{
        res.f <- "Unclassify"
    }

    return(res.f)
})
dat$tmp.group <- tmp.group
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/MSK123902_tumor_subgroup.pdf",useDingbats=F)
DimPlot(dat,group.by="tmp.group",cols=c("#ffc719","#2baf2b","#00a5dc","#ced7df"),reduction="tsne")
dev.off()




# 2021-11-23
dat$BMS.group <- "inte"
dat$BMS.group[which(dat$BMS.update<quantile(dat$BMS.update,0.25))] <- "low"
dat$BMS.group[which(dat$BMS.update>quantile(dat$BMS.update,0.75))] <- "high"
DefaultAssay(dat) <- "RNA"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")












###############################################################################################################################
# MSK data GSE123902 
# DTC Incipent and Macrometastasis

library(Seurat)
file <- dir("~/metastasis/data/verify/GSE123904/DTC2Macro/")
inte.list <- list()
for(i in 1:length(file)){
    tmp <- read.csv(paste0("~/metastasis/data/verify/GSE123904/DTC2Macro/",file[i]),stringsAsFactors=F,colClasses="character")

    if(i %in%c(3,4)){
        name <- strsplit(file[i],"_")[[1]][4]
    }else{
        name <- strsplit(file[i],"_")[[1]][3]
    }
    rownames(tmp) <- paste0("Cell_",tmp$X)
    tmp$X <- NULL
    tmp <- t(tmp)
    tmp.dat<- CreateSeuratObject(counts = tmp,  project = name ,min.cells = 3, min.features = 200)

    # MT gene 
    tmp.dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp.dat, pattern = "^MT.")

    # filter some cells
    dat = subset(x=tmp.dat,subset=nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
    dat = NormalizeData(object = dat)
    dat$group <- name
    dat$sample <- paste0(name,i)
    # add into inte.list 
    inte.list[[i]] <- dat
}

######################################################################################################################################
# integration samples 
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
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro.RDS")
# calculate sssGSEA for BMS and Brain and Lung 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","Lung_gene","Brain_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro_mod.RDS")



# check result 
################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro_mod.RDS")

dat$BMS.update <- mod$BMS_update_norm
dat$Lung.gene <- mod$Lung_gene_norm
dat$Brain.gene <- mod$Brain_gene_norm


# define cell sub group 
tmp.group <- apply(dat@meta.data[,c("Brain.gene","Lung.gene","BMS.update")],1,function(x){
    tmp <- c()
    res <- c()
    if(x[1]>quantile(dat$Brain.gene,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"group16")
    }
    if(x[2]>quantile(dat$Lung.gene,0.75)){
        tmp <- c(tmp,x[2])
        res <- c(res,"group09")
    }
    if(x[3]>quantile(dat$BMS.update,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"BMS.update")
    }

    if(length(res)>0){
        res.f <- res[which.max(tmp)]
    }else{
        res.f <- "Unclassify"
    }

    return(res.f)
})
dat$tmp.group <- tmp.group

tmp.res <- reshape2::melt(apply(table(dat$tmp.group,dat$group),1,function(x){x/sum(x)})[,c(1:3)])
colnames(tmp.res) <- c("group1","group2","percentage")

library(ggplot2)

cols <- c("#e9e8dd","#5ba4e5","#9fbb58")
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig4/MSKDTC2Macro.pdf",useDingbats=F)
ggplot(tmp.res,aes(x=group2,y=percentage,fill=group1,group=group1))+ geom_bar(stat="identity",position="dodge")+
    geom_text(aes(label=round(percentage,2), y=percentage+0.01), position=position_dodge(0.9), vjust=0)+
    scale_fill_manual(values=cols)+theme_bw()+ labs(y="percentage")
dev.off()



# 2021-11-23 
dat$BMS.group <- "inte"
dat$BMS.group[which(dat$BMS.update<quantile(dat$BMS.update,0.25))] <- "low"
dat$BMS.group[which(dat$BMS.update>quantile(dat$BMS.update,0.75))] <- "high"
DefaultAssay(dat) <- "RNA"
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123904_DTC2Macro.RDS")









# GSE14108 data analysis
# 2021-10-15
# estimate calculate purity
############################################################################################################################### 
library(estimate)
# calculate Brain metastasis sample purity 
#==============================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
# write.table(dat,file="GSE14108_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)
# filterCommonGenes(input.f="GSE14108_exp.txt",output.f="BM_10412genes.gct",id="GeneSymbol")
# estimateScore("BM_10412genes.gct", "/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GSE14108_estimate_score.gct", platform="affymetrix")

scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GSE14108_estimate_score.gct",skip = 2,header = T)


# ssGSEA 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Brain_gene","Lung_gene","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$purity <- as.numeric(as.vector(scores[4,-c(1,2)]))

mod$immune <- as.numeric(as.vector(scores[2,-c(1,2)]))


# pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_purity_lungsig.pdf",useDingbats=F)
# plot(mod$Lung_gene_norm,mod$purity,main="Lung like tumor with purity",xlab="Lung cancer signature",ylab="estimate purity score",pch=19)
# abline(lm(mod$purity~mod$Lung_gene_norm),col="red")
# # cor.tes()
# legend("topright",legend=paste0("rho= -",0.49," pvalue=",0.001))
# dev.off()


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_estimate_immune_lungsig.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,mod$immune,main="Lung like tumor with purity",xlab="Lung cancer signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$Lung_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.58," pvalue<",0.001))

plot(mod$Brain_gene_norm,mod$immune,main="Brain like tumor with purity",xlab="Brain cancer signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$Brain_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",-0.05," pvalue=",0.76))

plot(mod$BMS_update_norm,mod$immune,main="BMS with purity",xlab="BMS signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$BMS_update_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.47," pvalue=",0.01))

dev.off()





#================================================================================================================================
# 2021-10-15
# E-MTAB data analysi 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
# write.table(dat,file="E_MTAB_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)
# filterCommonGenes(input.f="E_MTAB_exp.txt",output.f="E_MTAB_genes.gct",id="GeneSymbol")
# estimateScore("E_MTAB_genes.gct", "/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/E_MTAB_estimate_score.gct", platform="affymetrix")

# calcualte 
scores <- read.table("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/E_MTAB_estimate_score.gct",skip = 2,header = T)
# ssGSEA 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Brain_gene","Lung_gene","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$purity <- as.numeric(as.vector(scores[4,-c(1,2)]))
mod$immune <- as.numeric(as.vector(scores[2,-c(1,2)]))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/E_MTAB_estimate_immune_lungsig.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,mod$immune,main="Lung like tumor with purity",xlab="Lung cancer signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$Lung_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.29," pvalue=",0.022))

plot(mod$Brain_gene_norm,mod$immune,main="Brain like tumor with purity",xlab="Brain cancer signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$Brain_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",-0.05," pvalue=",0.65))

plot(mod$BMS_update_norm,mod$immune,main="BMS with purity",xlab="BMS signature",ylab="estimate immune score",pch=19)
abline(lm(mod$immune~mod$BMS_update_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho= ",0.46," pvalue<",0.001))

dev.off()






#############################################################################################################
# 2021-8-18 
# to check gene in samples
#============================================================================================================
library(Seurat)

tumor.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
tumor.e@active.ident <- factor(tumor.e$tmp.group)

tumor.d <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
tumor.d@active.ident <- factor(tumor.d$tmp.group)

gse123902 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")

lcbm <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(lcbm) <- "RNA"
lcbm@active.ident <- factor(lcbm$group)

features <- c("CSF1","CCL20")

AverageExpression(tumor.e,assay="RNA",features=features)
AverageExpression(tumor.d,assay="RNA",features=features)
AverageExpression(gse123902,assay="RNA",features=features)
AverageExpression(lcbm,assay="RNA",features=features)






# 2021-9-9
# divide GSE123902 into 3 sample 
#============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
sub.dat <-  dat # subset(dat,cells=which(dat$sample=="LX681")) # LX255B  LX701 LX681
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
sub.dat$BMS.update <- mod[,4]
sub.dat$Brain.gene <- mod[,5]
sub.dat$Lung.gene <- mod[,6]

# define cell sub group 
tmp.group <- apply(sub.dat@meta.data[,c("Brain.gene","Lung.gene","BMS.update")],1,function(x){
    tmp <- c()
    res <- c()
    if(x[1]>quantile(sub.dat$Brain.gene,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"group16")
    }
    if(x[2]>quantile(sub.dat$Lung.gene,0.75)){
        tmp <- c(tmp,x[2])
        res <- c(res,"group09")
    }
    if(x[3]>quantile(sub.dat$BMS.update,0.75)){
        tmp <- c(tmp,x[1])
        res <- c(res,"BMS.update")
    }

    if(length(res)>0){
        res.f <- res[which.max(tmp)]
    }else{
        res.f <- "Unclassify"
    }

    return(res.f)
})
sub.dat$tmp.group <- tmp.group

saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")

sub.dat@active.ident <- factor(sub.dat$tmp.group)
AverageExpression(sub.dat,assays="RNA",features=c("CSF1","CCL2","SAA1"))


# saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor_LX681.RDS")




































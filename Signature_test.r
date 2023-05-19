#!/usr/bin/Rscript

# 2021-12-27
# BMS test program copy form Vesion5 
#========================================================================================================================================
# 1. GSE126548 test
# 2. single cell test
# 3. TCGA survival
# 4. cell line test 
# 5. TCGA satge test 
# 6. KEGG pathway analysis 
# 7. GSE68465 survival test
########################################################################################################################################
sigtest <- function(gene,name="gene_test"){

dir.create(paste0("/public/workspace/lily/Lung2Brain/Version6/Signature/",name))  # 2021-12-27 change this code 
respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Signature/",name,"/") # 2021-12-27 change this code 

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
library(survminer)
library(survival)
library(ggpubr)
library(dplyr)

# 0. make a mod
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,name,out=paste0(respath,name,".mod")) # make a mod file 

# 1. GSE126548
#=============================================================================================================
load('/public/workspace/lily/metastasis/data/verify/GSE126548/GSE126548_rpkm.RData')
gse126548_mod <- mod.analyze2(as.matrix(rpkm),name,respath,permN=0)
saveRDS(gse126548_mod,file=paste0(respath,'GSE126548_sig_mod.RDS'))

# plot result 
pdf(paste0(respath,'GSE126548_boxplot','.pdf'),useDingbats=F) # boxplot for GSE126548 (3 vs. 3)
boxplot(gse126548_mod[c(1:3),2],gse126548_mod[c(4:6),2],names=c('BM-','BM+'),main=paste0(name))
legend('topright',legend=paste('p =',round(wilcox.test(gse126548_mod[c(1:3),2],gse126548_mod[c(4:6),2])$p.value,5)),bty='n')
dev.off()



# 2. GSE123902 Single cell 
#==============================================================================================================
GSE123902 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")
GSE123902$group <- "Early"
GSE123902$group[which(GSE123902$sample=="LX675")] <- "Advanced"
gse123902_mod <- mod.analyze2(as.matrix(GSE123902[['RNA']]@data),name,respath,permN=0) # calculate mod
gse123902_mod <- data.frame(gse123902_mod)
saveRDS(gse123902_mod,file=paste0(respath,'GSE123902_sig_mod.RDS'))

# plot result 
gse123902_mod$group <- GSE123902$group
pdf(paste0(respath,'GSE123902_boxplot','.pdf'),useDingbats=F) # boxplot for GSE123902 (Single cell )
boxplot(gse123902_mod[,2]~group,data=gse123902_mod,FUN=median)
legend('topright',legend=paste('p =',round(wilcox.test(gse123902_mod[,2]~group,data=gse123902_mod,FUN=median)$p.value,5)),bty='n')
dev.off()






# 3. TCGA survival 
#===============================================================================================================
# use no normalized expression matrix (log2 trans)
# tmp <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_RNAseq_Exp.txt",header=T)
# rownames(tmp) <- tmp$sample
# tmp$sample <- NULL
# saveRDS(tmp,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
LUAD <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
LUAD_mod <- mod.analyze2(as.matrix(LUAD),name,respath,permN=0) # calculate mod
LUAD_mod <- data.frame(LUAD_mod)
saveRDS(LUAD_mod,file=paste0(respath,'TCGA_LUAD_sig_mod.RDS'))

# clinical information
# load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
# saveRDS(luad_ann,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RDS")
# rm(ann)
luad_ann <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RDS")

rownames(luad_ann) <- gsub("-",".",rownames(luad_ann))
luad_ann <- luad_ann[colnames(LUAD),]


luad_ann$type_five <- rep("M",nrow(luad_ann))
luad_ann[which(LUAD_mod[,2]>unname(quantile(LUAD_mod[,2],0.8))),"type_five"] <- "H"
luad_ann[which(LUAD_mod[,2]<unname(quantile(LUAD_mod[,2],0.2))),"type_five"] <- "L"
# LUAD half
luad_ann$type_two <- rep("L",nrow(luad_ann))
luad_ann[which(LUAD_mod[,2]>unname(quantile(LUAD_mod[,2],0.5))),"type_two"] <- "H"
# LUAD three quantile
luad_ann$type_three <- rep("M",nrow(luad_ann))
luad_ann[which(LUAD_mod[,2]>unname(quantile(LUAD_mod[,2],0.67))),"type_three"] <- "H"
luad_ann[which(LUAD_mod[,2]<unname(quantile(LUAD_mod[,2],0.33))),"type_three"] <- "L"

# plot result 
############ plot #################
pdf(paste0(respath,'sig_surv_LUAD','.pdf'),useDingbats=F)
#==================================
# plot OS time
# sig LUAD type : half 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survminer::surv_fit(surv~type_two,data=luad_ann)
p1 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_half",xlab=" Overall survival (months)")+labs(title="LUAD_sig_half")
# sig LUAD type : five 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survminer::surv_fit(surv~type_five,data=luad_ann)
p2 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_five",xlab=" Overall survival (months)")+labs(title="LUAD_sig_five")
# sig LUAD type : three 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survminer::surv_fit(surv~type_three,data=luad_ann)
p3 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" Overall survival (months)")+labs(title="LUAD_sig_three")

# plot RFS time
# sig1 LUAD type : half 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survminer::surv_fit(surv~type_two,data=luad_ann)
p4 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_half",xlab=" RF survival (months)")+labs(title="LUAD_sig_half")
# sig1 LUAD type1 : five 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survminer::surv_fit(surv~type_five,data=luad_ann)
p5 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_five",xlab=" RF survival (months)")+labs(title="LUAD_sig_five")
# sig1 LUAD type2 : three 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survminer::surv_fit(surv~type_three,data=luad_ann)
p6 <- ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" RF survival (months)")+labs(title="LUAD_sig_three")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

dev.off()


# 3.1 vasion and localvasion
#===============================================================================================================
invasion_mod <- mod.analyze2(as.matrix(LUAD),c("intraInvasion_cp","localInvasion_cp"),'/public/workspace/lily/MOD_file/',permN=0)
invasion_mod <- data.frame(invasion_mod)
saveRDS(invasion_mod,file=paste0(respath,'TCGA_LUAD_invasion_mod.RDS'))

pdf(paste0(respath,'TCGA_LUAD_vasion_plot','.pdf'),useDingbats=F) # dotplot
invasion_mod$sig_mod <- as.numeric(as.vector(LUAD_mod[,2]))
plot(invasion_mod$sig_mod,invasion_mod$intraInvasion_cp_norm)
legend('topright',legend=paste('rho =',round(cor.test(invasion_mod$sig_mod,invasion_mod$intraInvasion_cp_norm,method="spearman")$estimate,5),
", p=",round(cor.test(invasion_mod$sig_mod,invasion_mod$intraInvasion_cp_norm,method="spearman")$p.value,5)),bty='n')

plot(invasion_mod$sig_mod,invasion_mod$localInvasion_cp_norm)
legend('topright',legend=paste('rho =',round(cor.test(invasion_mod$sig_mod,invasion_mod$localInvasion_cp_norm,method="spearman")$estimate,5),
", p=",round(cor.test(invasion_mod$sig_mod,invasion_mod$localInvasion_cp_norm,method="spearman")$p.value,5)),bty='n')
dev.off()




# 4. cell line verfiy
# GSE14995
#=================================================================================================================
GSE14995 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14995/GSE14995_dat.RDS")
gse14995_mod <- mod.analyze2(as.matrix(GSE14995),name,respath,permN=0)
saveRDS(gse14995_mod,file=paste0(respath,'GSE14995_sig_mod.RDS'))

# plot result 
pdf(paste0(respath,'GSE14995_boxplot','.pdf'),useDingbats=F) # boxplot for GSE14995
boxplot(gse14995_mod[c(1:3),2],gse14995_mod[c(4:6),2],names=c('low invasion','high invasion'),main=paste0(name))
legend('topright',legend=paste('p =',round(wilcox.test(gse14995_mod[c(1:3),2],gse14995_mod[c(4:6),2])$p.value,5)),bty='n')
dev.off()




# 5. TCGA stage 
#================================================================================================================
# LUAD_mod
LUAD <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
# luad_mod
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
tcga.clin <- tcga.clin[which(rownames(tcga.clin)%in%colnames(LUAD)),]
all(colnames(LUAD)==rownames(tcga.clin))
LUAD <- LUAD[,rownames(tcga.clin)]
# filter some sample 
LUAD.f <- LUAD[,grep("Stage",tcga.clin$pathologic_stage)]
tcga.clin.f <- tcga.clin[which(rownames(tcga.clin)%in%colnames(LUAD.f)),]
#################################################################################
# Calculate BMS score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
luad.f.mod <- mod.analyze2(as.matrix(LUAD.f),name,respath,permN=0)
luad.f.mod <- data.frame(luad.f.mod)
luad.f.mod$stage <- "unknow"
luad.f.mod$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))] <- "early"
luad.f.mod$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))] <- "advance"
saveRDS(luad.f.mod,file=paste0(respath,'LUAD_fliter_sig_mod.RDS'))

# plot result 
pdf(paste0(respath,'LUAD_fliter_boxplot','.pdf'),useDingbats=F) # boxplot for TCGA LUAD filter 
boxplot(luad.f.mod[,2]~stage,data=luad.f.mod,FUN=median)
legend('topright',legend=paste('p =',round(wilcox.test(luad.f.mod[,2]~stage,data=luad.f.mod,FUN=median)$p.value,5)),bty='n')
dev.off()






# 6. KEGG pathway 
#====================================================================================================================
geneinfo <- bitr(gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
tmp1 <- ekk@result
tmp1 <- tmp1[which(tmp1$pvalue<0.05),]

# plot KEGG pathway

pdf(paste0(respath,'KEGG_pathway','.pdf'))
# sig 
tmp1$Description  <- factor(tmp1$Description,level=tmp1$Description[nrow(tmp1):1])
p_kegg <- ggplot(data=tmp1,aes(x=Description,y=-log(pvalue),fill=pvalue))+
	geom_bar(stat="identity",position='dodge',width=0.7)+ coord_flip()+
	theme_classic()+labs(x="Pathway",title='sig_KEGG')
print(p_kegg)
dev.off()





# 7. GSE68465 survival
#===================================================================================================================
GSE68465 <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
gse68465_info <- readRDS('/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS')

gse68465_mod <- mod.analyze2(as.matrix(GSE68465),name,respath,permN=0)
gse68465_mod <- as.data.frame(gse68465_mod)

res.gse68465 <- merge(gse68465_mod,gse68465_info,by="row.names")

# res.gse68465.f <- res.gse68465[-which(res.gse68465$Disease=="Normal"),]

res.ff <- res.gse68465
#================================================
# survival 
#================================================
library(survminer)
library(survival)

res.ff$type <- "MID"
res.ff$type[which(res.ff[,paste0(name,"_norm")]<unname(quantile(res.ff[,paste0(name,"_norm")],0.33)))] <- "Low"
res.ff$type[which(res.ff[,paste0(name,"_norm")]>unname(quantile(res.ff[,paste0(name,"_norm")],0.67)))] <- "High"
res.ff$OS_status <- 0 #change type info 
res.ff$OS_status[which(res.ff$Status=="Dead")] <- 1

table(res.ff$type)
res.ff$mths_to_last_clinical_assessment <- as.numeric(as.vector(res.ff$mths_to_last_clinical_assessment))
surv <- Surv(res.ff$mths_to_last_clinical_assessment,res.ff$OS_status)
km <- survfit(surv~type,data=res.ff)
pdf(paste0(respath,'GSE68465_survival','.pdf'))
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',
	title="GSE68465",xlab=" Overall survival (months)")+labs(title="GSE67465_sig_three")
dev.off()




}











# # annotation 
# # record how GSE123902 primary data come from 
# ###################################################################################################################################
# filelist <- grep("csv",dir("~/metastasis/data/verify/GSE123904/GSE123902"),value=T)

# inte.list <- list()
# for(i in 1:length(filelist)){
#     tmp <- as.data.frame(data.table::fread(paste0("~/metastasis/data/verify/GSE123904/GSE123902/",filelist[i])))
#     tmp$V1 <- paste0("Cell_",tmp$V1)
#     rownames(tmp) <- tmp$V1
#     tmp$V1 <- NULL
#     tmp.f <- t(tmp)

#     # create a seurat object 
#     tmp.obj <- CreateSeuratObject(counts=tmp.f,project = strsplit(filelist[i],"_")[[1]][4])
#     tmp.obj$sample <- as.vector(strsplit(filelist[i],"_")[[1]][3])
#     tmp.obj$group <- as.vector(strsplit(filelist[i],"_")[[1]][4])
#     tmp.obj = NormalizeData(object = tmp.obj)
# 	tmp.obj <- FindVariableFeatures(object = tmp.obj)

#     inte.list[i] <- tmp.obj
# }

# integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
# inte <- IntegrateData(anchorset = integration.anchors)
# #FindVariableFeatures
# inte <- FindVariableFeatures(inte)
# ##Scaling the integrateda
# all.genes <- rownames(inte)
# inte <- ScaleData(inte, features = all.genes)
# #PCA
# inte <- RunPCA(inte)
# #cluster
# inte <- FindNeighbors(inte)
# inte <- FindClusters(inte)
# #TSNE
# # if Umap can not use
# inte <- RunTSNE(inte)
# saveRDS(inte,file="/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")

# # FeaturePlot(inte,features=c("EPCAM","EGFR"),label=T)
# # show cluster 6 maybe tumor cell beacuse of EPCAM expression 
# # https://www.nature.com/articles/s41591-019-0750-6/figures/7  this show primary sample 
# # and this data do not have therapy data, do not have other metastasis and 2021-6-23 lly copy it to Version5/MSK_scRNA/
# #=======================================================================================================================================
# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")
# sub.dat <- subset(dat,cells=which(dat$seurat_clusters==6&dat$group=="PRIMARY"))


























































































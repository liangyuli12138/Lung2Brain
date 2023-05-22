#!/usr/bin/Rscript

# 2021-4-30
# this program is used to plot Fig2
#========================================================================================
# 0. use sig_test.r change
# 

source("/public/workspace/lily/Lung2Brain/Version5/Fig2/sig_test.r")
sigtest(name="BMS_update")






#########################################################################################
# other some validate 

# 1. GSE8894 
#==================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/ann.RDS")
ann.f <- ann[which(ann$Type=="Adenocarcinoma"),]
dat.f <- dat[,which(colnames(dat)%in%rownames(ann.f))]
source('~/software/ssGSEA/ssgseaMOD.r')

mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
mod$Recurrent <- ann.f$Recurrent

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/GSE8894_recurrent.pdf",useDingbats=F)
boxplot(BMS_update_norm~Recurrent,data=mod,FUN=median,outline=F,ylim=c(0,1))
legend("topright",legend=paste0("P.value=",wilcox.test(BMS_update_norm~Recurrent,data=mod)$p.value))
dev.off()




# 2. TCGA LUAD early stage verify 
#=========================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage","RFS.time","RFS","OS.time","OS",
"new_neoplasm_event_type","new_tumor_event_after_initial_treatment","days_to_new_tumor_event_after_initial_treatment")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
##########################################################################################################################
tcga.f <- dat[,which(colnames(dat)%in%rownames(tcga.clin)[which(tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))])]

# calculate BMS score 
###########################################################################################################################
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

sample.high <- rownames(mod)[which(mod$BMS_update_norm>=quantile(mod$BMS_update_norm,0.67))]
sample.low <- rownames(mod)[which(mod$BMS_update_norm<=quantile(mod$BMS_update_norm,0.33))]

tcga.clin.f <- tcga.clin[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)

# OS time
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/TCGA_early_stage_surv.pdf",useDingbats=F)
surv <- Surv(tcga.clin.f$OS.time,tcga.clin.f$OS)
km <- survfit(surv~group,data=tcga.clin.f)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),
surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Overall survival (months)")+labs(title="Early stage LUAD")

# RFS time 
surv <- Surv(tcga.clin.f$RFS.time,tcga.clin.f$RFS)
km <- survfit(surv~group,data=tcga.clin.f)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),
surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Progression free survival (months)")+labs(title="Early stage LUAD")
dev.off()





# 2021-8-11
# 2021-8-12
# try to use BMS to calculate Distant Metastasis and non-Metastasis sample 
#============================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage","RFS.time","RFS","OS.time","OS",
"new_neoplasm_event_type","new_tumor_event_after_initial_treatment","days_to_new_tumor_event_after_initial_treatment")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
# get sample names
early <- rownames(tcga.clin)[which(tcga.clin$new_tumor_event_after_initial_treatment=="NO"&tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))]
earlyM <- rownames(tcga.clin)[which(tcga.clin$new_neoplasm_event_type=="Distant Metastasis"&tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))]
#  get TCGA data and clinical data 
tcga.f <- dat[,which(colnames(dat)%in%c(early,earlyM))]
tcga.clin.sub <- tcga.clin[which(rownames(tcga.clin)%in%c(early,earlyM)),]


# get sample names
# nMet <- rownames(tcga.clin)[which(tcga.clin$new_tumor_event_after_initial_treatment=="NO")]
# Met <- rownames(tcga.clin)[which(tcga.clin$new_neoplasm_event_type=="Distant Metastasis")]
# #  get TCGA data and clinical data 
# tcga.f <- dat[,which(colnames(dat)%in%c(nMet,Met))]
# tcga.clin.sub <- tcga.clin[which(rownames(tcga.clin)%in%c(nMet,Met)),]

# tranlation Metastasis info 
tcga.clin.sub$Met <- 0
tcga.clin.sub$Met[which(tcga.clin.sub$new_tumor_event_after_initial_treatment=="YES")] <- 1

#===================================================================================================================================
# calculate BMS and plot survival
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

# mod.rs <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
# mod <- mod.rs[which(rownames(mod.rs)%in%rownames(mod)),]
# all(rownames(mod)==colnames(tcga.f))

sample.high <- rownames(mod)[which(mod$BMS_update_norm>=quantile(mod$BMS_update_norm,0.75))]
sample.low <- rownames(mod)[which(mod$BMS_update_norm<=quantile(mod$BMS_update_norm,0.25))]

tcga.clin.f <- tcga.clin.sub[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/early_TCGA_Metastasis_surv_25p.pdf",useDingbats=F)
surv <- Surv(tcga.clin.f$RFS.time,tcga.clin.f$Met)
km <- survfit(surv~group,data=tcga.clin.f)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),
surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Progression free survival (months)")+labs(title="Early stage LUAD")
dev.off()


# try to use boxplot to check result 
# however not significant
tcga_clin <- tcga.clin.sub[which(rownames(tcga.clin.sub)%in%rownames(mod)),]
tcga_clin <- tcga_clin[rownames(mod),]
all(rownames(tcga_clin)==rownames(mod))

mod$group <- tcga_clin$Met











# 3. other survival validate  
# GSE68465 
# 2021-4-30 seem this just have trend but not signifcant
#=========================================================================================================================

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
info <- readRDS('/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS')

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
res.ff <- merge(mod,info,by="row.names")

#================================================
# survival 
#================================================
library(survminer)
library(survival)

res.ff$type <- "MID"
res.ff$type[which(res.ff$BMS_update_norm<unname(quantile(res.ff$BMS_update_norm,0.33)))] <- "Low"
res.ff$type[which(res.ff$BMS_update_norm>unname(quantile(res.ff$BMS_update_norm,0.67)))] <- "High"
res.ff$OS_status <- 0 #change type info 
res.ff$OS_status[which(res.ff$Status=="Dead")] <- 1

res.ff <- res.ff[-which(res.ff$type=="MID"),]
res.ff$OS.time <- as.numeric(as.vector(res.ff$months_to_last_contact_or_death))
surv <- Surv(res.ff$OS.time,res.ff$OS_status)
km <- survfit(surv~type,data=res.ff)
pairwise_survdiff(Surv(res.ff$OS.time,res.ff$OS_status)~type,data=res.ff)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,surv.median.line='hv')$plot+scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))

dev.off()






# 2021-9-26
# add CCL20 to analysis
library(survminer)
library(survival)

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
info <- readRDS('/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS')
survinfo <- info[colnames(dat),]

survinfo$CCL20 <- as.numeric(dat["CCL20",])
survinfo$type <- "MID"
survinfo$type[which(survinfo$CCL20<unname(quantile(survinfo$CCL20,0.5)))] <- "Low"
survinfo$type[which(survinfo$CCL20>unname(quantile(survinfo$CCL20,0.5)))] <- "High"
survinfo$OS_status <- 0 #change type info 
survinfo$OS_status[which(survinfo$Status=="Dead")] <- 1

info.f <- survinfo[-which(survinfo$type=="MID"),]
info.f$OS.time <- as.numeric(as.vector(info.f$months_to_last_contact_or_death))
surv <- Surv(info.f$OS.time,info.f$OS_status)
km <- survfit(surv~type,data=info.f)
pairwise_survdiff(Surv(info.f$OS.time,info.f$OS_status)~type,data=info.f)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,surv.median.line='hv')$plot+scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))






























# 4. make sure this signature can be find in LCBM group
# 2021-5-10
#============================================================================================================
library(Seurat)
recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.5)

    tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10)
	return(tmp_dat)
}

tumor.d <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
tumor.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
tumor.d <- recluster(tumor.d)
tumor.e <- recluster(tumor.e)

# calculate BMS 
source('~/software/ssGSEA/ssgseaMOD.r')
mod.d <- mod.analyze2(as.matrix(tumor.d[["RNA"]]@data),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod.d <- data.frame(mod.d)

mod.e <- mod.analyze2(as.matrix(tumor.e[["RNA"]]@data),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod.e <- data.frame(mod.e)

# 
tumor.d$BMS_update <- mod.d[,2]
tumor.e$BMS_update <- mod.e[,2]


pdf("./tmp.pdf")
DimPlot(tumor.d,label=T,reduction="tsne")
FeaturePlot(tumor.d,features="BMS_update",reduction="tsne",min.cutoff=0.5)
DimPlot(tumor.e,label=T,reduction="tsne")
FeaturePlot(tumor.e,features="BMS_update",reduction="tsne",min.cutoff=0.5)

DimPlot(tumor.d,label=T)
FeaturePlot(tumor.d,features="BMS_update",min.cutoff=0.6)
DimPlot(tumor.e,label=T)
FeaturePlot(tumor.e,features="BMS_update",min.cutoff=0.6)
dev.off()









# 5. why choose C4 and C7 
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Signature/LCBM_hallmark.RDS")
# mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Signature/LCBM_metastasis_pathway.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Signature/pur_LCBM_clust.RDS")

# mod.f <- mod[,22:42]
mod.f <- mod[,51:100]
mod.f$Cluster <- paste0("C",dat$seurat_clusters)
tmp.res <- aggregate(.~Cluster,data=mod.f,FUN=median)
rownames(tmp.res) <- tmp.res$Cluster
tmp.res$Cluster <- NULL
res.f <- t(tmp.res)

table(apply(res.f,1,function(x){colnames(res.f)[order(x,decreasing=T)]})[1:2,])
# calculate hallmark and top2 cluster 
# > table(apply(res.f,1,function(x){colnames(res.f)[order(x,decreasing=T)]})[1:2,])

#  C0  C1 C10  C2  C4  C6  C7  C8  C9
#  16   1  16   1  29   1  21   7   8





# add a asian data 
# GSE31210
# use shell to get sample info 
library(GEOquery)
gset <- getGEO( "GSE31210", getGPL = F )
Gset <- gset[[1]]
pdata<-pData(Gset)
sampleinfo <- pdata[,c(2,11,62,64,52,53,54,57)]
colnames(sampleinfo) <- c("GEO","Age","Stage","RFS","OS.time","RFS.time","OS","Gene_fusion")
sampleinfo$Age <- as.numeric(gsub("age \\(years\\): |age: ","",sampleinfo$Age))
sampleinfo$RFS[which(sampleinfo$RFS=="not relapsed")] <- 0
sampleinfo$RFS[which(sampleinfo$RFS=="relapsed")] <- 1
sampleinfo$OS[which(sampleinfo$OS=="alive")] <- 0
sampleinfo$OS[which(sampleinfo$OS=="dead")] <- 1
info.f <- sampleinfo[-c(227:246),]
saveRDS(info.f,file="~/metastasis/data/verify/GSE31210/sampleinfo.RDS")

# data 
dat <- read.table("~/metastasis/data/verify/GSE31210/GSE31210_series_matrix.txt",sep="\t",header=T,comment.char="!")
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/GPL570.txt",sep="\t",header=T)
tmp.dat <- merge(tmp,dat,by.x="ID",by.y="ID_REF")
tmp.dat.f <- tmp.dat[-grep("//",tmp.dat$Gene.Symbol),-1]
tmp.res <- aggregate(.~Gene.Symbol,data=tmp.dat.f,FUN=mean)
res.f <- tmp.res[-c(1:25),]
rownames(res.f) <- res.f$Gene.Symbol
res.f$Gene.Symbol <- NULL

saveRDS(res.f,file="~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")


###################################################################################################################################
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
info <- readRDS("~/metastasis/data/verify/GSE31210/sampleinfo.RDS")
# dat.f <- dat[,which(colnames(dat)%in%rownames(info))]
# dat.f <- dat.f[,rownames(info)]
# saveRDS(dat.f,file="~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
all(rownames(mod)==rownames(info))
info$BMS <- mod[,2]
# make survival 
library(survminer)
library(survival)

info$type <- "MID" # a bit question is MID is the best
info$type[which(info$BMS<unname(quantile(info$BMS,0.33)))] <- "Low"
info$type[which(info$BMS>unname(quantile(info$BMS,0.67)))] <- "High"
# change type
info$OS.time <- as.numeric(as.vector(info$OS.time))
info$RFS.time <- as.numeric(as.vector(info$RFS.time))
info$OS <- as.numeric(as.vector(info$OS))
info$RFS <- as.numeric(as.vector(info$RFS))


info <- info[-which(info$type=="MID"),]
surv <- Surv(info$OS.time,info$OS)
km <- survfit(surv~type,data=info)
pairwise_survdiff(Surv(info$OS.time,info$OS)~type,data=info)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/GSE31210_surv.pdf",useDingbats=F)
ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Overall survival (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")

surv <- Surv(info$RFS.time,info$RFS)
km <- survfit(surv~type,data=info)
pairwise_survdiff(Surv(info$RFS.time,info$RFS)~type,data=info)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" RFS (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")

dev.off()








#================================================================================================================================================
# 2021-9-26
# data analysis for CCL20 
#================================================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
info <- readRDS("~/metastasis/data/verify/GSE31210/sampleinfo.RDS")

# make survival 
library(survminer)
library(survival)

info$CCL20 <- as.numeric(dat["CCL20",])
info$type <- "MID" # a bit question is MID is the best
info$type[which(info$CCL20<unname(quantile(info$CCL20,0.33)))] <- "Low"
info$type[which(info$CCL20>unname(quantile(info$CCL20,0.67)))] <- "High"
# change type
info$OS.time <- as.numeric(as.vector(info$OS.time))
info$RFS.time <- as.numeric(as.vector(info$RFS.time))
info$OS <- as.numeric(as.vector(info$OS))
info$RFS <- as.numeric(as.vector(info$RFS))

info.f <- info[-which(info$type=="MID"),]
surv <- Surv(info.f$OS.time,info.f$OS)
km <- survfit(surv~type,data=info.f)
pairwise_survdiff(Surv(info.f$OS.time,info.f$OS)~type,data=info.f)

ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Overall survival (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")


surv <- Surv(info.f$RFS.time,info.f$RFS)
km <- survfit(surv~type,data=info.f)
pairwise_survdiff(Surv(info.f$RFS.time,info.f$RFS)~type,data=info.f)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" RFS (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")





















###########################################################################################################################
# COX analysis for TCGA LUAD 
# 2021-5-17
#==========================================================================================================================
luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","EGFR","KRAS","MET","ALK_translocation","ROS1_translocation","ERBB2","BRAF",
	"age_at_initial_pathologic_diagnosis","OS.time","OS","RFS.time","RFS","gender","pathologic_stage","vital_status",
	"Cnncl_mt_n_KRAS_EGFR_ALK_RET_ROS1_BRAF_ERBB2_HRAS_NRAS_AKT1_MAP2")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
tcga.clin.f <- merge(tcga.clin,luad_mod,by="row.names")

# BMS 
# tcga.clin.f$BMS_type <- "MID"
# tcga.clin.f$BMS_type[which(tcga.clin.f$BMS_test_norm>quantile(tcga.clin.f$BMS_test_norm,0.67))] <- "High"
# tcga.clin.f$BMS_type[which(tcga.clin.f$BMS_test_norm<quantile(tcga.clin.f$BMS_test_norm,0.33))] <- "Low"


# age
tcga.clin.f$age_at_initial_pathologic_diagnosis -> tcga.clin.f$Age
tcga.clin.f$age_at_initial_pathologic_diagnosis <- NULL

# mutation use 0 1 2 to subtitued
tcga.clin.f$EGFR_mut <- 1
tcga.clin.f$EGFR_mut[which(tcga.clin.f$EGFR=="none")] <- 0
tcga.clin.f$EGFR_mut[which(tcga.clin.f$EGFR=="")] <- NA


# mutation use 0 1 2 to subtitued
tcga.clin.f$KRAS_mut <- 1
tcga.clin.f$KRAS_mut[which(tcga.clin.f$KRAS=="none")] <- 0
tcga.clin.f$KRAS_mut[which(tcga.clin.f$KRAS=="")] <- NA


# mutation use 0 1 2 to subtitued
tcga.clin.f$MET_mut <- 1
tcga.clin.f$MET_mut[which(tcga.clin.f$MET=="none")] <- 0
tcga.clin.f$MET_mut[which(tcga.clin.f$MET=="")] <- NA

tcga.clin.f$ALK_trans <- 1
tcga.clin.f$ALK_trans[which(tcga.clin.f$ALK_translocation=="none")] <- 0
tcga.clin.f$ALK_trans[which(tcga.clin.f$ALK_translocation=="")] <- NA

# ROS1 
tcga.clin.f$ROS1_trans <- 1
tcga.clin.f$ROS1_trans[which(tcga.clin.f$ROS1_translocation=="none")] <- 0
tcga.clin.f$ROS1_trans[which(tcga.clin.f$ROS1_translocation=="")] <- NA

# BRAF
tcga.clin.f$BRAF_mut <- 1
tcga.clin.f$BRAF_mut[which(tcga.clin.f$BRAF=="none")] <- 0
tcga.clin.f$BRAF_mut[which(tcga.clin.f$BRAF=="")] <- NA

# ERBB2 
tcga.clin.f$ERBB2_mut <- 1
tcga.clin.f$ERBB2_mut[which(tcga.clin.f$ERBB2=="none")] <- 0
tcga.clin.f$ERBB2_mut[which(tcga.clin.f$ERBB2=="")] <- NA

# gender use 0 1 to change 
tcga.clin.f$sex <- 0
tcga.clin.f$sex[which(tcga.clin.f$gender=="MALE")] <- 1


# Stage use 0  1 2 3 4 to change 

tcga.clin.f$stage <- NA
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB","Stage IIIA","Stage IIIB"))] <- "Low stage"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IV"))] <- "High stage"
tcga.clin.f$stage <- factor(tcga.clin.f$stage,levels=c("Low stage","High stage"))


library(ggplot2)
library(survival)
library(survminer)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/TCGA_Coxph.pdf",useDingbats=F)

model <- coxph(Surv(OS.time, OS) ~  BMS_update_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+BRAF_mut+ERBB2_mut+stage, data = tcga.clin.f)
ggforest(model, data = tcga.clin.f)

dev.off()











#####################################################################################################################################################
# 2021-5-17
# TCGA LUAD 
# early stage and late stage 
#====================================================================================================================================================
luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","EGFR","KRAS","MET","ALK_translocation","ROS1_translocation","ERBB2","BRAF",
	"age_at_initial_pathologic_diagnosis","OS.time","OS","RFS.time","RFS","gender","pathologic_stage")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
tcga.clin.f <- merge(tcga.clin,luad_mod,by="row.names")

##########################################################################################################################
tcga.clin.f$stage <- NA
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))] <- "Early stage"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))] <- "Advanced stage"
res.f <- tcga.clin.f[-which(is.na(tcga.clin.f$stage)),]
# plot result 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/TCGA_LUAD_early_late.pdf",useDingbats=F)
ggplot(dat=res.f ,aes(x=stage,y=BMS_update_norm,fill=stage))+stat_boxplot(geom = "errorbar",width=0.15)+geom_boxplot(notch = T,outlier.shape = NA) + theme_classic()+
    annotate("text",x=2,y=1,label=wilcox.test(BMS_update_norm~stage,data=res.f,FUN=median)$p.value,parse=T,color="blue")
dev.off()






# GSE123902
# 2021-5-18
#====================================================================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
GSE123902 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")
GSE123902$group <- "Early"
GSE123902$group[which(GSE123902$sample=="LX675")] <- "Advanced"
gse123902_mod <- mod.analyze2(as.matrix(GSE123902[['RNA']]@data),"BMS_update","/public/workspace/lily/MOD_file/",permN=0) # calculate mod
gse123902_mod <- data.frame(gse123902_mod)
gse123902_mod$group <- GSE123902$group

library(ggridges)
library(ggplot2)

# boxplot


# TCGA extravasion 
#===================================================================================================================================================
tmp1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig2/TCGA_LUAD_invasion_mod.RDS")
tmp2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
res <- data.frame(BMS=tmp2$BMS_update_norm,localvasion=tmp1$localInvasion_cp_norm,extra=tmp1$intraInvasion_cp_norm)

# plot result 
library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/BMS_extravasion_cor.pdf",useDingbats=F)
p<-ggplot(res,aes(x=BMS,y=extra)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Extravasion score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "spearman",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
p1 <- ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
print(p1)
dev.off()


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig2/BMS_localvasion_cor.pdf",useDingbats=F)
p<-ggplot(res,aes(x=BMS,y=localvasion)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("localvasion score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "spearman",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()

































# 2021-9-9 
# calculate in GSE131907 data by sample 
#===================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/useTumor.RDS")
tmp.res <- list()
# calculate by each sample 
samplelist <- as.vector(unique(dat$Sample))
for(i in 1:length(samplelist)){
    tmp.sub <- subset(dat,cells=which(dat$Sample==samplelist[i]))
    DefaultAssay(tmp.sub) <- "RNA"
    # calculate BMS score 
    source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
    mod <- mod.analyze2(as.matrix(tmp.sub[['RNA']]@data),c("BMS_update","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
    mod <- as.data.frame(mod)
    tmp.res[[i]] <- mod
    names(tmp.res)[i] <- samplelist[i]
}

saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_Sig.RDS")




# tmp
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_Sig.RDS")
res <- sapply(tmp.res,function(x){median(x[,6])})


















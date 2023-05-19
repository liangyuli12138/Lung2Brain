
# 2022-1-3 
# this program is used to analysis BMS and verify BMS 
# 0. plot score sketch
# 1. analysis BMS score and paired data 
# 2. early stage TCGA survival
# 3. early stage MFS 
# 4. recurrent sample verify 
# 5. GSE31210 survival analysis
# 6. coxph analysis 
# 7. replot extravasion plot 
# 2022-5-26 add GSE200563 analysis 
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================
library(Seurat)
library(monocle)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
# get Paired tumor cells 
subdat <- subset(dat,cells=which(dat$orig.ident%in%c("Pair_BM","Pair_Lung")&dat$celltype.refine=="Tumor"))
# now re-inte to find genes 
inte.list <- list() 
samplelist <- unique(subdat$orig.ident) 

for(i in 1:length(samplelist)){
	tmp <- subset(subdat,cells=which(subdat$orig.ident==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
inte <- IntegrateData(anchorset = integration.anchors)
#FindVariableFeatures
inte <- FindVariableFeatures(inte)
var.gene <- VariableFeatures(inte) # get gene

# also calulate BMS score
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(inte[['RNA']]@data),c("BMS_V6","HPSC_C5"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Paired_LCBM_tumor_BMS_mod.RDS")
inte$BMS_score <- mod$BMS_V6_norm

# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/Paired_LCBM_tumor_BMS_mod.RDS")

#==========================================================================================================================================
# trans into monocle cds data 
data <- as(as.matrix(inte@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = inte@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- var.gene[1:1500]
cds <- setOrderingFilter(cds, ordering_genes)
# pseudotime 
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)

saveRDS(cds,file="/public/workspace/lily/Lung2Brain/Version6/Data/Pair_Tumor_1.5K_monocle.RDS")
plot_cell_trajectory(cds,color_by="BMS_score")
plot_cell_trajectory(cds,color_by="type_group")



#==========================================================================================================================================
#==========================================================================================================================================
# 2. do some analysis about TCGA early stage sample and MFS 
#==========================================================================================================================================
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
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_V6"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

sample.high <- rownames(mod)[which(mod$BMS_V6_norm>=quantile(mod$BMS_V6_norm,0.67))]
sample.low <- rownames(mod)[which(mod$BMS_V6_norm<=quantile(mod$BMS_V6_norm,0.33))]

tcga.clin.f <- tcga.clin[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)

# OS time
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/TCGA_early_stage_surv.pdf",useDingbats=F)
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



###########################################################################################################################
# MFS 
# 2022-1-7
###########################################################################################################################
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

# tranlation Metastasis info 
tcga.clin.sub$Met <- 0
tcga.clin.sub$Met[which(tcga.clin.sub$new_tumor_event_after_initial_treatment=="YES")] <- 1

#===================================================================================================================================
# calculate BMS and plot survival
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_V6"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

# annotation info 
sample.high <- rownames(mod)[which(mod$BMS_V6_norm>=quantile(mod$BMS_V6_norm,0.67))]
sample.low <- rownames(mod)[which(mod$BMS_V6_norm<=quantile(mod$BMS_V6_norm,0.33))]

tcga.clin.f <- tcga.clin.sub[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/early_TCGA_Metastasis_surv_33p.pdf",useDingbats=F)
surv <- Surv(tcga.clin.f$RFS.time,tcga.clin.f$Met)
km <- survfit(surv~group,data=tcga.clin.f)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),
surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Progression free survival (months)")+labs(title="Early stage LUAD")
dev.off()










#==========================================================================================================================================
#==========================================================================================================================================
# 4. do some analysis recurrent sample data 
# GSE8894
#==========================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/ann.RDS")
ann.f <- ann[which(ann$Type=="Adenocarcinoma"),]
dat.f <- dat[,which(colnames(dat)%in%rownames(ann.f))]
source('~/software/ssGSEA/ssgseaMOD.r')

mod <- mod.analyze2(as.matrix(dat.f),c("BMS_V6"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
mod$Recurrent <- ann.f$Recurrent

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/GSE8894_recurrent.pdf",useDingbats=F)
boxplot(BMS_V6_norm~Recurrent,data=mod,FUN=median,outline=F,ylim=c(0,1))
legend("topright",legend=paste0("P.value=",wilcox.test(BMS_V6_norm~Recurrent,data=mod)$p.value))
dev.off()







###################################################################################################################################
# GSE31210
# asian survival
# early stage 
#==================================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
info <- readRDS("~/metastasis/data/verify/GSE31210/sampleinfo.RDS")
# dat.f <- dat[,which(colnames(dat)%in%rownames(info))]
# dat.f <- dat.f[,rownames(info)]
# saveRDS(dat.f,file="~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_V6"),"/public/workspace/lily/MOD_file/",permN=0)
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

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/GSE31210_surv.pdf",useDingbats=F)
ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" Overall survival (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")

surv <- Surv(info$RFS.time,info$RFS)
km <- survfit(surv~type,data=info)
pairwise_survdiff(Surv(info$RFS.time,info$RFS)~type,data=info)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,xlim=c(0,3000),xscale='d_m',break.time.by =304.37,xlab=" RFS (months)")$plot+
scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))+labs(title="Early stage LUAD")

dev.off()









###########################################################################################################################
# COX analysis for TCGA LUAD 
# 2022-1-7
#==========================================================================================================================
# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("BMS_V6"),'/public/workspace/lily/MOD_file/',permN=0)
# mod <- data.frame(mod)
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")

luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
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
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))] <- "Low stage"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))] <- "High stage"
tcga.clin.f$stage <- factor(tcga.clin.f$stage,levels=c("Low stage","High stage"))


library(ggplot2)
library(survival)
library(survminer)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/TCGA_Coxph.pdf",useDingbats=F)

model <- coxph(Surv(OS.time, OS) ~  BMS_V6_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+BRAF_mut+ERBB2_mut+stage, data = tcga.clin.f)
ggforest(model, data = tcga.clin.f)

dev.off()


#================================================================================================================================================
# single 
# forest plot 

library(ggplot2)
library(survival)
library(survminer)
covariates <- colnames(tcga.clin.f)[c(19,21:30)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(RFS.time, RFS)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = tcga.clin.f)})

fdat <- data.frame(factors = names(sapply(univ_models,function(x){summary(x)$conf.int[3]})),
    exp = round(as.numeric(sapply(univ_models,function(x){summary(x)$conf.int[1]})),2),
    lower = round(as.numeric(sapply(univ_models,function(x){summary(x)$conf.int[3]})),2),
    upper =  round(as.numeric(sapply(univ_models,function(x){summary(x)$conf.int[4]})),2),
    exp_conf = NA,
    pvalue = round(as.numeric(sapply(univ_models,function(x){summary(x)$logtest[3]})),4),
    stringsAsFactors=F)
fdat$exp_conf <- paste0(fdat$exp," [",fdat$lower,",",fdat$upper,"]")

# model <- coxph(Surv(RFS.time, RFS) ~  BMS_V6_norm+ERBB2_mut+stage, data = tcga.clin.f)
# ggforest(model, data = tcga.clin.f)

library(forestplot)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/TCGA_RFS_single_Coxph.pdf",useDingbats=F)
forestplot(labeltext=fdat, graph.pos=3,
    mean=c(fdat$exp),
    lower=c(fdat$lower), upper=c(fdat$upper),
    zero=1,
    cex=0.9, lineheight = "auto",
    colgap=unit(8,"mm"),
    #箱子大小，线的宽度
    lwd.ci=2, boxsize=0.2,
    #箱线图两端添加小竖线，高度
    ci.vertices=TRUE, ci.vertices.height = 0.25)
dev.off()

tcga.clin.f$stage <- NA
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))] <- "I"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in% c("Stage II","Stage IIA","Stage IIB"))] <- "II"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB"))] <- "III"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IV"))] <- "IV"





###########################################################################################################################
# replot extravasion plot for TCGA 
# 2022-1-7
#==========================================================================================================================
# LUAD <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# invasion_mod <- mod.analyze2(as.matrix(LUAD),c("intraInvasion_cp","localInvasion_cp"),'/public/workspace/lily/MOD_file/',permN=0)
# invasion_mod <- data.frame(invasion_mod)
# saveRDS(invasion_mod,file=paste0('/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_invasion_mod.RDS'))

# TCGA extravasion 
#==========================================================================================================================
tmp1 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_invasion_mod.RDS")
tmp2 <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
res <- data.frame(BMS=tmp2$BMS_V6_norm,localvasion=tmp1$localInvasion_cp_norm,extra=tmp1$intraInvasion_cp_norm)

# plot result 
library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/BMS_extravasion_cor.pdf",useDingbats=F)
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


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/BMS_localvasion_cor.pdf",useDingbats=F)
p<-ggplot(res,aes(x=BMS,y=localvasion)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("localvasion score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "spearman",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()










#===========================================================================================================================
# 2022-5-26
# GSE200563 
# plot result for BMS in this data

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
tmp.res <- readRDS("~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")

mod.gse <- as.data.frame(mod.analyze2(as.matrix(tmp.res$GSE_Data),"BMS_V6","/public/workspace/lily/MOD_file/",permN=0))
mod.tcga <- as.data.frame(mod.analyze2(as.matrix(tmp.res$TCGA_Data),"BMS_V6","/public/workspace/lily/MOD_file/",permN=0))


plotdat <- data.frame(group=c(rep("TCGA",nrow(mod.tcga)),as.vector(tmp.res$GSE_info$group)),BMS=c(mod.tcga[,2],mod.gse[,2]))
# make a data frame to plot violin 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/TCGA_GSE200563_BMS.pdf",useDingbats=F)
boxplot(BMS~group,plotdat,outline=F,ylim=c(0,1))
dev.off()







#===========================================================================================================================
# 2022-6-7
# GSE200563 
# do correlation for BMS and time metastasis to brain

tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("MLUNG")),]
dat <- tmp.dat[,rownames(sampleinfo)]

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_V6"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)



mod$braintime <- sampleinfo$metastasis_intervals_to_brain

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/GSE200563_BMS_braintime_correlation.pdf",useDingbats=F)
plot(mod$BMS_V6_norm,as.numeric(mod$braintime),main=" GSE200563 (MLUNG,N=27)")
abline(lm(as.numeric(mod$braintime)~mod$BMS_V6_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$BMS_V6_norm,as.numeric(mod$braintime),method="spearman")$estimate,
    " P =",round(cor.test(mod$BMS_V6_norm,as.numeric(mod$braintime),method="spearman")$p.value,2))
)
dev.off()








#=====================================================================================================================================================
# plot skect 
# 2022-11-15
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/All_tumorcells.pdf",useDingbats=F)
DimPlot(dat,cells.highlight=colnames(dat)[which(dat$celltype.refine=="Tumor")],reduction="tsne",raster=F)
dev.off()


LCBM <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Signature/LCBM_tumor_pur.RDS")
nMLUAD <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Signature/nMLUAD_tumor_pur.RDS")

#subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(6,8,10,14)))
# subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(2,8,6)))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig2/LCBM_nMLUAD_subset.pdf",useDingbats=F)
DimPlot(LCBM,cells.highlight=colnames(LCBM)[which(LCBM$seurat_clusters%in%c(6,8,10,14))],reduction="tsne",raster=F)
DimPlot(nMLUAD,cells.highlight=colnames(nMLUAD)[which(nMLUAD$seurat_clusters%in%c(2,8,6))],reduction="tsne",raster=F)
dev.off()






































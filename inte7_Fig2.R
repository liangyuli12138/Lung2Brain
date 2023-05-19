

#=======================================================================================================
# 2020-10-6
# coxph plot
#=======================================================================================================
library(ggplot2)
library(survival)
library(survminer)










#=======================================================================================================
# sketch fig to know how to find gene marker 
#=======================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
BM_t <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
LC_t <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/sketchplot.pdf",useDingbats=F)
DimPlot(dat,cells.highlight=colnames(dat)[which(dat$maliganant=="tumor")])
DimPlot(BM_t,cells.highlight=colnames(BM_t)[which(BM_t$seurat_clusters%in%c(6,9))])
DimPlot(LC_t,cells.highlight=colnames(LC_t)[which(LC_t$seurat_clusters%in%c(0,11))])
dev.off()


















#=======================================================================================================
# TCGA samples 

load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","EGFR","KRAS","MET","ALK_translocation","ROS1_translocation","ERBB2","BRAF",
	"age_at_initial_pathologic_diagnosis","OS.time","OS","RFS.time","RFS","gender","pathologic_stage","vital_status")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
tcga.clin.f <- merge(tcga.clin,luad_mod,by="row.names")

#================================
# BMS type classify
tcga.clin.f$BMS_type <- "MID"
tcga.clin.f$BMS_type[which(tcga.clin.f$BMS_test_norm>quantile(tcga.clin.f$BMS_test_norm,0.67))] <- "High"
tcga.clin.f$BMS_type[which(tcga.clin.f$BMS_test_norm<quantile(tcga.clin.f$BMS_test_norm,0.33))] <- "Low"


#================================
# age
tcga.clin.f$age_at_initial_pathologic_diagnosis -> tcga.clin.f$Age
tcga.clin.f$age_at_initial_pathologic_diagnosis <- NULL





#=================================================================================
# model <- coxph( Surv(OS.time, OS) ~  Age +gender + BMS_type + Canonical_mut_in_KRAS_EGFR_ALK , data = tcga.clin.f )
# model <- coxph( Surv(OS.time, OS) ~  BMS_type , data = tcga.clin.f )
# model
# ggforest(model, data = tcga.clin.f)


#================================================================================
# maybe have a bit question 
#
#================================================================================
# BMS use BMS norm to subtitued
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
tcga.clin.f$sex[which(tcga.clin.f$gender=="FEMALE")] <- 1


# Stage use 0  1 2 3 4 to change 

tcga.clin.f$stage <- NA
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))] <- "I"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage II","Stage IIA","Stage IIB"))] <- "II"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB"))] <- "III"
tcga.clin.f$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IV"))] <- "IV"


model <- coxph( Surv(OS.time, OS) ~  BMS_test_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+stage , data = tcga.clin.f )

model <- coxph( Surv(RFS.time, RFS) ~  BMS_test_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+stage , data = tcga.clin.f )


# more factor 
model <- coxph( Surv(OS.time, OS) ~  BMS_test_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+ERBB2_mut+BRAF_mut+stage , data = tcga.clin.f )

model <- coxph( Surv(RFS.time, RFS) ~  BMS_test_norm+Age+sex+EGFR_mut+KRAS_mut+MET_mut+stage , data = tcga.clin.f )


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/cox.pdf",useDingbats=F)
ggforest(model, data = tcga.clin.f)
dev.off()




#========================================================================================================
# modify TCGA STAGE info 
# maybe some III patients should add into IV 
# do not know how to find 
#========================================================================================================
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]













#===============================================================================================
# KEGG pathway anlysis 
#===============================================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig
gene <- readRDS("/public/workspace/lily/Lung2Brain/inte7/BMS_gene.RDS")
geneinfo <- bitr(gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
tmp1 <- ekk@result
tmp1 <- tmp1[which(tmp1$pvalue<0.05),]

# plot KEGG pathway

pdf(paste0("/public/workspace/lily/Lung2Brain/inte7/Fig/",'BMS_KEGG','.pdf'))
# sig 
tmp1$Description  <- factor(tmp1$Description,level=tmp1$Description[nrow(tmp1):1])
ggplot(data=tmp1,aes(x=Description,y=-log(pvalue),fill=pvalue))+
	geom_bar(stat="identity",position='dodge',width=0.7)+ coord_flip()+
	theme_classic()+labs(x="Pathway",title='sig_KEGG')
dev.off()




#======================================================================================
# 2020-12-17
# GSE131907 use BM to verify
#======================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat<- subset(dat,cells=which(dat$Sample_Origin%in%c("mBrain","nLung","tL/B","tLung")&dat$Cell_type%in%c("Epithelial cells")))
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)





#====================================================================================
# GSE149246
# use Cell line data to verify
#####################################################################################
for OBJ in `ls ./*.txt`;
do
cut -f 1,6,7 $OBJ > ${OBJ}.tmp
done
# paste ./*.tmp > counts.txt
# do some prepare and make final count.txt 

tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE149246/count.txt",sep="\t",header=T)
anno <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",sep="\t",header=T)
tmp.res <- merge(anno[,c(1,2)],tmp.dat,by.x="Gene.stable.ID",by.y="Geneid")
# transform into expression data 
################################################################################################################
# tmp.res[,2] <- NULL
# 
tmp <- t(apply(tmp.res[,3:11],1,function(x){x/x[1]}))*10^3
res.f <- apply(tmp,2,function(x){x/sum(x)})*10^6
res.f <- data.frame(res.f)
res.f$gene_name <-tmp.res$Gene.name
res.data <- aggregate(.~gene_name,data=res.f,FUN=median)
rownames(res.data) <- res.data$gene_name
res.data[,c(1,2)] <- NULL
saveRDS(as.matrix(res.data),file="/public/workspace/lily/metastasis/data/verify/GSE149246/GSE149246_expr.RDS")


#==============================================================================================================
# ssGSEA calculate 
# just use sample 1 to 6 ,because sample 7 and sample 8 is shHSF1
###############################################################################################################
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE149246/GSE149246_expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[,1:6]),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE149246_boxplot.pdf")
boxplot(mod[1:2,2],mod[3:6,2])
legend("topright",legend=paste0("P=",round(wilcox.test(mod[1:2,2],mod[3:6,2])$p.value,digits=3)))
dev.off()





###################################################################################################################
# local vasion and extravasion calculate 
#==================================================================================================================
# TCGA data verify
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('BMS_test',"intraInvasion_cp","localInvasion_cp"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)


# just plot intravasion res
###################################################################################################################

library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/BMS_extravasion_cor.pdf",useDingbats=F)
p<-ggplot(mod,aes(x=BMS_test_norm,y=intraInvasion_cp_norm)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Extravasion score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()



#################################################################################################################
# tumor and adjcent normal 
# GSE 40419
########################################
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE40419/GSE40419_LC-87_RPKM_expression.txt",sep="\t",header=T)
anno <- read.delim("/public/workspace/lily/REF/INDEX-hg19/anno/RefseqID_symbol.txt",sep="\t",header=T)
tmp.res <- merge(anno[,c(1,2)],tmp.dat,by.x="RefSeq.IDs",by.y="accession")
tmp.res[,c(1,3:7)] <- NULL
res.f <- aggregate(.~Approved.symbol,data=tmp.res,FUN=median)
rownames(res.f) <- res.f$Approved.symbol
res.f$Approved.symbol <- NULL
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/GSE40419/GSE40419expr.RDS")

# ssGSEA calculate 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE40419/GSE40419expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
mod$type <- c(rep("LC",87),rep("nor",77))



##########################################
# density plot use ggplot2 
#=========================================
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE40419_tumor_normal_density.pdf")
ggplot(data=mod,aes(x=BMS_test_norm,fill=type)) + geom_density(alpha=0.4) + theme_classic()
dev.off()



pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE40419_tumor_normal_boxplot.pdf")
boxplot(BMS_test_norm~type,data=mod,FUN=median)
legend("topright",legend=paste0("P=",round(wilcox.test(BMS_test_norm~type,data=mod,FUN=median)$p.value,digits=3)))
dev.off()





########################################################################################################################
# DisGENT signature verify
# ###### Fri Jan 15 17:14:32 CST 2021
#=======================================================================================================================
library(Seurat)
BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
subset1 <- subset(BM,cells=which(BM$seurat_clusters%in%c(6,9)))
subset2 <- subset(LC,cells=which(LC$seurat_clusters%in%c(0,11)))
tmp <- merge(subset1,subset2)
# tmp$type_group <- "LC"
# tmp$type_group[which(tmp$orig.ident%in%c("A20190305","A20190312","T-Bsc1"))] <- "LCBM"
tmp@active.ident <- as.factor(tmp$type_group)
geneset_tmp <- FindMarkers(tmp,ident.1="LCBM",ident.2="LC",assay="RNA",logfc.threshold = 0.1)
geneset_tmp$pct <- geneset_tmp$pct.1-geneset_tmp$pct.2
geneset_tmp$pct1 <- geneset_tmp$pct.1/geneset_tmp$pct.2


tmp.res <- geneset_tmp[which(geneset_tmp$avg_logFC>1&geneset_tmp$p_val_adj<0.01&geneset_tmp$pct>0.2&geneset_tmp$pct1>2),]
tmp.f <- tmp.res[-grep("MT-",rownames(tmp.res))]
tmp.f <- tmp.f[order(tmp.f$avg_logFC,decreasing=T),]
write.table(tmp.f[,c("avg_logFC"),drop=F],file="./BMS_rank.txt",row.names=T,col.names=T,quote=F)
# then add "Gene" 









#=======================================================================================================================
# GSE68465
# surv 
#=======================================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
info <- readRDS('/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS')

mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)
mod <- as.data.frame(mod)

res.ff <- merge(mod,info,by="row.names")

#res.ff <- res.f[-which(res.f$Disease=="Normal"),]
#================================================
# survival 
#================================================
library(survminer)
library(survival)

res.ff$type <- "MID"
res.ff$type[which(res.ff$BMS_test_norm<unname(quantile(res.ff$BMS_test_norm,0.33)))] <- "Low"
res.ff$type[which(res.ff$BMS_test_norm>unname(quantile(res.ff$BMS_test_norm,0.67)))] <- "High"
res.ff$OS_status <- 0 #change type info 
res.ff$OS_status[which(res.ff$Status=="Dead")] <- 1

table(res.ff$type)
res.ff$mths_to_last_clinical_assessment <- as.numeric(as.vector(res.ff$mths_to_last_clinical_assessment))
surv <- Surv(res.ff$mths_to_last_clinical_assessment,res.ff$OS_status)
km <- survfit(surv~type,data=res.ff)
pairwise_survdiff(Surv(res.ff$mths_to_last_clinical_assessment,res.ff$OS_status)~type,data=res.ff)
pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,surv.median.line='hv')$plot+scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))
dev.off()















#===================================================================================================================
# use local vasion and extravasion data to verify Sangsung data 
#2021-3-11
#===================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/useTumor.RDS")
gene <- list(c('TGFB1','CSF1','EGF','VEGFA','PTGS2','EREG','ANGPT2','MMP1','MMP2','MMP3','MMP10'))
dat <- AddModuleScore(dat,features=gene,name="extra")

# add module score 
boxplot(extra1~Sample_Origin,data=dat@meta.data,FUN=median)


# use ssgsea 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c('intraInvasion_cp',"localInvasion_cp"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
mod$type  <- dat$Sample_Origin
mod$sample <- dat$Samples
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/GSE131907/useTumor_extravasion_mod.RDS")









#==========================================================================================================================
# 2021-3-11
# use TCGA LUAD data and GSE14108 to verify BMS
# 2021-4-15
# just use TCGA sample to verfiy  early and late stage 
#==========================================================================================================================

load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # TCGA LUAD Data
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
##########################################################################################################################
dat.f <- apply(dat,2,function(x){x/x[8176]}) # use ACTB to adjust 
tcga.f <- dat.f[,which(colnames(dat.f)%in%tcga.clin$sampleID[which(tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))])]

# ggplot2 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_LUAD_BM.pdf",useDingbats=F)
ggplot(dat=mod,aes(x=type,y=BMS_test_norm,fill=type))+ geom_violin()+geom_boxplot(width=.1) + theme_classic()
dev.off()

###########################################################################################################################
# try to use none normalized TCGA LUAD samples
# not good result 
#==========================================================================================================================
# dat <- read.delim2("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_RNAseq_Exp.txt",sep="\t",header=T)
# rownames(dat) <- dat$sample
# dat$sample <- NULL
# tmp <- apply(dat,2,function(x){as.numeric(as.vector(x))})
# rownames(tmp) <- rownames(dat)
# dat <- tmp 

# tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
# tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]
# rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
# ##########################################################################################################################
# dat.f <- apply(dat,2,function(x){x/x[8176]}) # use ACTB to adjust 
# tcga.f <- dat.f[,which(colnames(dat.f)%in%rownames(tcga.clin)[which(tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))])]

# # Lung cancer brain metastsis samples 
# gse <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS") # GSE data 
# gse.f <- apply(gse,2,function(x){x/x[149]})


# res.f <- merge(tcga.f,gse.f,by="row.names")
# rownames(res.f) <- res.f$Row.names
# res.f$Row.names <- NULL









##################################################################################################################################
# 2021-4-13
# check TCGA stage I difference 
#=================================================================================================================================
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # TCGA LUAD Data

tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage","RFS.time","RFS","OS.time","OS")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
##########################################################################################################################
tcga.f <- dat[,which(colnames(dat)%in%tcga.clin$sampleID[which(tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))])]

# calculate BMS score 
#=========================================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

sample.high <- rownames(mod)[which(mod$BMS_test_norm>=quantile(mod$BMS_test_norm,0.67))]
sample.low <- rownames(mod)[which(mod$BMS_test_norm<=quantile(mod$BMS_test_norm,0.33))]

tcga.clin.f <- tcga.clin[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)

surv <- Surv(tcga.clin.f$RFS.time,tcga.clin.f$RFS)
km <- survfit(surv~group,data=tcga.clin.f)
# pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,xlim=c(0,1095),break.x.by=365)$plot+
    scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))
    
    
    
    +
    xlim(c(0,1100))
# dev.off()





#================================================================================================================================
# GSE8894 test recurrent and BMS score 
#################################################################################################################################
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/ann.RDS")
ann.f <- ann[which(ann$Type=="Adenocarcinoma"),]
dat.f <- dat[,which(colnames(dat)%in%rownames(ann.f))]
source('~/software/ssGSEA/ssgseaMOD.r')

mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
mod$Recurrent <- ann.f$Recurrent

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE8894_recurrent.pdf",useDingbats=F)
boxplot(BMS_test_norm~Recurrent,data=mod,FUN=median,outline=F,ylim=c(0,1))
legend("topright",legend=paste0("P.value=",wilcox.test(BMS_test_norm~Recurrent,data=mod)$p.value))
dev.off()













# 2021-4-14
# deal with GSE123902 
###################################################################################################################################
filelist <- grep("csv",dir("~/metastasis/data/verify/GSE123904/GSE123902"),value=T)

inte.list <- list()
for(i in 1:length(filelist)){
    tmp <- as.data.frame(data.table::fread(paste0("~/metastasis/data/verify/GSE123904/GSE123902/",filelist[i])))
    tmp$V1 <- paste0("Cell_",tmp$V1)
    rownames(tmp) <- tmp$V1
    tmp$V1 <- NULL
    tmp.f <- t(tmp)

    # create a seurat object 
    tmp.obj <- CreateSeuratObject(counts=tmp.f,project = strsplit(filelist[i],"_")[[1]][4])
    tmp.obj$sample <- as.vector(strsplit(filelist[i],"_")[[1]][3])
    tmp.obj$group <- as.vector(strsplit(filelist[i],"_")[[1]][4])
    tmp.obj = NormalizeData(object = tmp.obj)
	tmp.obj <- FindVariableFeatures(object = tmp.obj)

    inte.list[i] <- tmp.obj
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
inte <- FindClusters(inte)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
saveRDS(inte,file="/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")

# FeaturePlot(inte,features=c("EPCAM","EGFR"),label=T)
# show cluster 6 maybe tumor cell beacuse of EPCAM expression 
# https://www.nature.com/articles/s41591-019-0750-6/figures/7  this show primary sample 
#=======================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")
sub.dat <- subset(dat,cells=which(dat$seurat_clusters==6&dat$group=="PRIMARY"))

# sub.dat <- AddModuleScore(sub.dat,features=gene,name="BMS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
mod$group <- "Early"
mod$group[which(sub.dat$sample=="LX675")] <- "Advanced"

mod <- mod[order(mod$group,mod$BMS_test_norm),]
mod$cellnum <- c( c(1:length(which(mod$group=="Advanced")))/length(which(mod$group=="Advanced")),c(1:length(which(mod$group=="Early")))/length(which(mod$group=="Early")) )


# point plot
# library(ggplot2)
# ggplot(data=mod,aes(y=BMS_test_norm,x=cellnum,col=group))+geom_point(aes(size=BMS_test_norm))+
# 	#scale_x_continuous(breaks=seq(1, max(mod$cellnum), 1000)) +
# 	theme_classic()
pdf('/public/workspace/lily/Lung2Brain/inte7/Fig/GSE123902_LungT_advanced.pdf',useDingbats=F,width=20)
library(ggridges)
library(ggplot2)
ggplot(mod) +
  geom_density_ridges_gradient(
    aes(x = BMS_test_norm, y = group, fill = group,height = ..density..), 
    scale = 1, rel_min_height = 0.01) +
  theme_ridges(font_size = 13, grid = TRUE) 

dev.off()

















################################################################################################################################################
# add TCGA lung early and late 
# 2021-4-16
#===============================================================================================================================================
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # TCGA LUAD Data ## dat
colnames(dat) <- gsub("-",".",colnames(dat))
# luad_mod
tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
tcga.clin <- tcga.clin[which(rownames(tcga.clin)%in%colnames(dat)),]
all(colnames(dat)==rownames(tcga.clin))
# filter some sample 
dat.f <- dat[,grep("Stage",tcga.clin$pathologic_stage)]
tcga.clin.f <- tcga.clin[which(rownames(tcga.clin)%in%colnames(dat.f)),]
###############################################################################################################################################
# Calculate BMS score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)

mod$stage <- "unknow"
mod$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))] <- "early"
mod$stage[which(tcga.clin.f$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))] <- "advance"


# do not use unknow sample 
# plot result by violin 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_LUAD_early_late.pdf",useDingbats=F)
ggplot(dat=mod,aes(x=stage,y=BMS_test_norm,fill=stage))+geom_boxplot(notch = T,outlier.shape = NA) + theme_classic()+
    annotate("text",x=2,y=1,label=wilcox.test(BMS_test_norm~stage,data=mod,FUN=median)$p.value,parse=T,color="blue")

dev.off()

# boxplot(BMS_test_norm~stage,data=mod,FUN=median)










#===========================================================================================================================================
# TCGA stage I and stage II sample 
# to check result 
#===========================================================================================================================================
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # TCGA LUAD Data

tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage","RFS.time","RFS","OS.time","OS")]
rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
##########################################################################################################################
tcga.f <- dat[,which(colnames(dat)%in%tcga.clin$sampleID[which(tcga.clin$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))])]

# calculate BMS score 
#=========================================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tcga.f),c("BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)

sample.high <- rownames(mod)[which(mod$BMS_test_norm>=quantile(mod$BMS_test_norm,0.67))]
sample.low <- rownames(mod)[which(mod$BMS_test_norm<=quantile(mod$BMS_test_norm,0.33))]


tcga.clin.f <- tcga.clin[c(sample.high,sample.low),]
tcga.clin.f$group <- "High"
tcga.clin.f$group[which(rownames(tcga.clin.f)%in%sample.low)] <- "Low"

library(survminer)
library(survival)



pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_early_stage_surv.pdf",useDingbats=F)
# OS time 
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
    





















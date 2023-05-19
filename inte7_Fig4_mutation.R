
#===================================================================================================
# 2020-12-17
# Mutation to be a new Fig 
#===================================================================================================
# 




#========================================================================================================================================
# another analysis about KRAS mutation 
#========================================================================================================================================
# first calculate 
#==========================================
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
a <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_broad_LUAD',sep='\t',header=T)
a$sample <- gsub('-','.',a$sample)
mut <- a[which(a$sample%in%rownames(luad_mod)),]
mod <- luad_mod[which(rownames(luad_mod)%in%mut$sample),]

gene <- unique(mut$gene)
res <- data.frame(gene)
p_value <- c()
mut_mean <- c()
no_mut_mean <- c()
percent <- c()
for(i in 1:length(gene)){
	as.vector(unique(mut[which(mut$gene==gene[i]),]$sample)) -> samp
	per <- length(unique(mut[which(mut$gene==gene[i]),]$sample))/length(unique(mut$sample))
	which(rownames(mod)%in%as.vector(samp)) -> y
	p <- wilcox.test(mod[y,2],mod[-y,2])$p.value
	y_mean <- mean(mod[y,2])
	n_mean <- mean(mod[-y,2])
	p_value <- c(p_value,p/2)
	mut_mean <- c(mut_mean,y_mean)
	no_mut_mean<- c(no_mut_mean,n_mean)
	percent <- c(percent,per)
}
cbind(res,mut_mean,no_mut_mean,percent,p_value) -> res
res$p.adj <- p.adjust(res$p_value)
res-> luad_res1

#=================================================================================================
# calculate rank
#=======================================
res.f <- luad_res1[which(luad_res1$p_value<0.01),]
res.ff <- res.f[which(res.f$percent>0.2),]
res.ff$varation <- (res.ff$mut_mean-res.ff$no_mut_mean)
res.ff <- res.ff[order(res.ff$varation,decreasing=T),]
# rank by P and percentage 

# res.ff
# res.f <- res.f[order(res.f$p_value,res.f$percent,decreasing=c(T,T)),]

#================================================================================================
# barplot to show difference 
#================================================================================================
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_LUAD_mutation_gene.pdf")
barplot((res.ff$varation-min(res.ff$varation))/(max(res.ff$varation)-min(res.ff$varation)),
    names=as.vector(res.ff$gene),
    ylab="normalized BMS score varation")
dev.off()


#=================================================================================================
# calculate CCLE data to verify 
#=================================================================================================
# TP53
dat <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct",skip=2,sep="\t",header=T)
dat$Name <- NULL
dat.f <- aggregate(.~Description,data=dat,FUN=median)
saveRDS(dat.f,file="/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS") # save expression data 

# calculate TP53
sampleinfo <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_TP53_mut.txt",sep="\t",header=T)
expr <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
# filter cell line 
expr.f <- expr[,which(colnames(expr)%in%sampleinfo$CCLE_ID)]

# calculate ssGSEA result 
#===============================================================================================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(expr.f),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=1000)
mod <- as.data.frame(mod)
mod$CCLE_ID <- rownames(mod)
res <- merge(mod,sampleinfo,by="CCLE_ID")




#===============================================================================================
# KRAS
# calculate TP53
sampleinfo <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_TP53_KRAS_mut.txt",sep="\t",header=T)
expr <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
# filter cell line 
expr.f <- expr[,which(colnames(expr)%in%sampleinfo$CCLE_ID)]

# calculate ssGSEA result 
#===============================================================================================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(expr.f),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=1000)
mod <- as.data.frame(mod)
mod$CCLE_ID <- rownames(mod)
res <- merge(mod,sampleinfo,by="CCLE_ID")



#===============================================================================================
# ready to cox plot 
#
#===============================================================================================
# make sampleinfo a matrix 
#===============================================================================================
sampleinfo <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_all_mut.txt",sep="\t",header=T)
sampleinfo$TP53 <- 0
sampleinfo$KRAS <- 0
sampleinfo$STK11 <- 0
sampleinfo$KEAP1 <- 0

for( i in 1:nrow(sampleinfo) ){
    if(length(grep("TP53",sampleinfo$Mutation[i]))>0){sampleinfo[i,"TP53"] <- 1}
    if(length(grep("KRAS",sampleinfo$Mutation[i]))>0){sampleinfo[i,"KRAS"] <- 1}
    if(length(grep("STK11",sampleinfo$Mutation[i]))>0){sampleinfo[i,"STK11"] <- 1}
    if(length(grep("KEAP1",sampleinfo$Mutation[i]))>0){sampleinfo[i,"KEAP1"] <- 1}
}

expr <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
# filter cell line 
expr.f <- expr[,which(colnames(expr)%in%sampleinfo$CCLE_ID)]
#===============================================================================================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(expr.f),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=1000)
mod <- as.data.frame(mod)
mod$CCLE_ID <- rownames(mod)
res <- merge(mod,sampleinfo,by="CCLE_ID")



#==============================================================================================
# set some result for cox analysis 
#==============================================================================================
res.f <- res[,c("CCLE_ID","BMS_test_norm","TP53","KRAS","KEAP1","STK11")]
res.f$sensors <- 1
write.table(res.f,file="/public/workspace/lily/metastasis/data/verify/CCLE/cox_sample_info.txt",sep="\t",row.names=F,col.names=T,quote=F)
write.table(expr.f,file="/public/workspace/lily/metastasis/data/verify/CCLE/cox_expr.txt",sep="\t",row.names=T,col.names=T,quote=F)

















#============================================================================================================================
# get all Cell line sample mutation info 
#
#============================================================================================================================
dat <- read.table("~/metastasis/data/verify/CCLE/CCLE_Mut_info.txt",sep=" ",header=T)
Cellline <- unique(as.vector(dat$Broad_ID))
res <- matrix(NA,ncol=2)
for(i in 1:length(Cellline)){
    mut_gene <- paste(as.vector(dat$Hugo_Symbol[which(dat$Broad_ID==Cellline[i])]),collapse=",")
    res <- rbind(res,c(Cellline[i],mut_gene))
}

#============================================================================================================================
# use LUAD sample info 
# run in local R studio 
#============================================================================================================================

setwd("F:\\lly\\Lung2Brain\\CCLE")
sampleinfo <- read.table("./CCLE_all_mut.txt",sep="\t",header = T)
options(stringsAsFactors = F)
samples <- unique(as.vector(sampleinfo$CCLE_ID))
sampleinfo$TP53 <- 1
sampleinfo$KRAS <- 1
sampleinfo$RP1L1 <- 1
sampleinfo$PCLO <- 1
sampleinfo$TTN <- 1
sampleinfo$MUC16 <- 1
sampleinfo$CSMD3 <- 1
sampleinfo$LRP1B <- 1

# tp53
tmp <- read.delim2("./TP53.txt",sep="\t")
sampleinfo$TP53[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# kras
tmp <- read.delim2("./KRAS.txt",sep="\t")
sampleinfo$KRAS[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# rp1l1
tmp <- read.delim2("./RP1L1.txt",sep="\t")
sampleinfo$RP1L1[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# pclo
tmp <- read.delim2("./PCLO.txt",sep="\t")
sampleinfo$PCLO[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# ttn 
tmp <- read.delim2("./TTN.txt",sep="\t")
sampleinfo$TTN[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# muc16
tmp <- read.delim2("./MUC16.txt",sep="\t")
sampleinfo$MUC16[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# csmd3
tmp <- read.delim2("./CSMD3.txt",sep="\t")
sampleinfo$CSMD3[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0

# lrp1b
tmp <- read.delim2("./LRP1B.txt",sep="\t")
sampleinfo$LRP1B[-which(as.vector(sampleinfo$CCLE_ID)%in%tmp$Tumor.Sample.Barcode)] <- 0


res.f <- sampleinfo
res.f$Mutation <- NULL
write.table(res.f,file="./8mutation.txt",sep = "\t",row.names = T)




#=================================================================================================================================
# calculate final 
# 2020-11-27
#=================================================================================================================================
sampleinfo <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/8mutation.txt",sep="\t",header=T)
expr <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
# filter cell line 
expr.f <- expr[,which(colnames(expr)%in%sampleinfo$CCLE_ID)]

source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(expr.f),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=1000)
mod <- as.data.frame(mod)
mod$CCLE_ID <- rownames(mod)
res <- merge(mod,sampleinfo,by="CCLE_ID")

#==================================================
# use BIT in 202.195.187.5
# run cox analysis 
#==================================================
res.f <- res[,-c(2,4)]
res.f$cencors <- 1
write.table(res.f,file="/public/workspace/lily/metastasis/data/verify/CCLE/8mut_cox_info.txt",sep="\t",quote=F,row.names=F)

library(survival)
library(survminer)
cox.res <- coxph(Surv(BMS_test_norm, cencors) ~ TP53+KRAS+RP1L1+PCLO+TTN+MUC16+CSMD3+LRP1B, data = res.f)
ggforest(cox.res)
# res.f$TP53 <- factor(res.f$TP53,level=c(1,0))
# res.f$KRAS <- factor(res.f$KRAS,level=c(1,0))
# res.f$RP1L1 <- factor(res.f$RP1L1,level=c(1,0))
# res.f$PCLO <- factor(res.f$PCLO,level=c(1,0))
# res.f$TTN <- factor(res.f$TTN,level=c(1,0))
# res.f$MUC16 <- factor(res.f$MUC16,level=c(1,0))
# res.f$CSMD3 <- factor(res.f$CSMD3,level=c(1,0))
# res.f$LRP1B <- factor(res.f$LRP1B,level=c(1,0))

#=========================================================================================
# and do single factor analysis 
#=========================================================================================
surv <- Surv(res.f$BMS_test_norm,res.f$cencors)
km <- survfit(surv~TP53,data=res.f)
ggsurvplot(km,pvalue=T)
cox.TP53 <- coxph(Surv(BMS_test_norm, cencors) ~ TP53, data = res.f)






#=========================================================================================
# use glm to plot 
# logsitic 
# the result should use forestplot 
# 2020-12-10
#=========================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/8mut_cox_info.txt",sep="\t",header=T)

library(survival)
library(survminer)

tmp <- glm(BMS_test_norm~TP53+KRAS+RP1L1+PCLO+TTN+MUC16+CSMD3+LRP1B,data=dat)
tmp.CI <- confint(tmp,level=0.9)  # should install MASS packages 

#========================================================================================
# make a input data 
#========================================================================================
tmp.res <- as.data.frame(summary(tmp)$coefficients[,c(1,4)])
colnames(tmp.res) <- c("value","Pvalue")
tmp.res <- cbind(tmp.res,tmp.CI)
colnames(tmp.res) <- c("value","Pvalue","5%CI","95%CI")

#========================================================================================
library(forestplot)
# plot 
final.dat <- tmp.res
final.dat$OR <- paste0(final.dat$value," (",final.dat$`5%CI`," ",final.dat$`95%CI`,")")
final.dat$gene <- rownames(final.dat) 
final.dat  <- final.dat[-1,]
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/forest_Gene_mut.pdf",useDingbats=F)
forestplot(rownames(final.dat),exp(final.dat$value),exp(final.dat$`5%CI`),exp(final.dat$`95%CI`), zero = 1, cex  = 2,
           lineheight = "auto",
           xlab = "Lab axis txt")
dev.off()

forestplot(as.matrix(final.dat[,c("value")]),final.dat$`5%CI`,final.dat$`95%CI`)







#===============================================================================================
# 2020-12-18
# calculate maf plot result 
#===============================================================================================
#########KRAS##########
#library(g3viz)
library(ggplot2)
library(maftools)
#########LUAD-C1##########
load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData')
luad_mod <- luad_mod[order(luad_mod[,2],decreasing=T),]
a <- luad_mod[,2]
nM <- luad_mod[which(luad_mod[,2]<unname(quantile(a,0.33))),]
M <- luad_mod[which(luad_mod[,2]>unname(quantile(a,0.67))),]
rownames(nM) <- gsub("\\.",'-',rownames(nM))
rownames(M) <- gsub("\\.","-",rownames(M))

maffile <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_LUAD.maf',sep='\t',header=T)
maffile[which(maffile$Tumor_Sample_Barcode%in%rownames(nM)),]-> maf_nM
maffile[which(maffile$Tumor_Sample_Barcode%in%rownames(M)),]-> maf_M
maf_nM$Tumor_Sample_Barcode <- factor(maf_nM$Tumor_Sample_Barcode)
maf_M$Tumor_Sample_Barcode <- factor(maf_M$Tumor_Sample_Barcode)

maf_n <-  read.maf(maf_nM)
maf_y <- read.maf(maf_M)


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/KRAS_maftools.pdf",useDingbats=F)
lollipopPlot(maf_n,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
printCount=T,showMutationRate = TRUE,showDomainLabel=F)

lollipopPlot(maf_y,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,printCount=T,
showMutationRate = TRUE,showDomainLabel=F) 
dev.off()




pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TP53_maftools.pdf",useDingbats=F)
lollipopPlot(maf_n,gene = "TP53",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
printCount=T,showMutationRate = TRUE,showDomainLabel=F)

lollipopPlot(maf_y,gene = "TP53",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,printCount=T,repel=T,
showMutationRate = TRUE,showDomainLabel=F) 
dev.off()





#==========================================================================================================
# 2020-12-18
# verify if KRAS mutation TCGA samples have KRAS pathway activate
# 
#==========================================================================================================
maffile <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_LUAD.maf',sep='\t',header=T)
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData")
load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData')
KRAS_mut_sample <- unique(as.vector(maffile$Tumor_Sample_Barcode[which(maffile$Hugo_Symbol=="KRAS")]))

source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- as.data.frame(mod)
mod$KRAS_type <- "non"
mod$KRAS_type[which(rownames(mod)%in%gsub("-",".",KRAS_mut_sample))] <- "mut"
aggregate(c(HALLMARK_KRAS_SIGNALING_UP_norm-HALLMARK_KRAS_SIGNALING_DN_norm)~KRAS_type,data=mod,FUN=median)

# 2021-2-19
# use RAS pathway activation signature from BMC article
#===========================================================================================================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("RAS.pathway.act"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- as.data.frame(mod)
mod$KRAS_type <- "non"
mod$KRAS_type[which(rownames(mod)%in%gsub("-",".",KRAS_mut_sample))] <- "mut"
aggregate(RAS.pathway.act_norm~KRAS_type,data=mod,FUN=median)

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/RAS_pathway_KRAS_mut.pdf",useDingbats=F)
boxplot(RAS.pathway.act_norm~KRAS_type,data=mod,FUN=median,outline=F,ylim=c(0,1))
legend("topright",paste0("P = ",wilcox.test(RAS.pathway.act_norm~KRAS_type,data=mod,FUN=median)$p.value))
dev.off()



# pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Hallmark_KRAS_activate.pdf",useDingbats=F)
# boxplot(c(HALLMARK_KRAS_SIGNALING_UP_norm-HALLMARK_KRAS_SIGNALING_DN_norm)~KRAS_type,data=mod,outline=F)
# legend("topright",paste0("P = ",wilcox.test(c(HALLMARK_KRAS_SIGNALING_UP_norm-HALLMARK_KRAS_SIGNALING_DN_norm)~KRAS_type,data=mod)$p.value))
# dev.off()

# use BMS mod 
#==========================================================================================================
# mod$BMS_type <- "median"
# mod$BMS_type[which(luad_mod[,2]<unname(quantile(luad_mod[,2],0.33)))] <- "Low"
# mod$BMS_type[which(luad_mod[,2]>unname(quantile(luad_mod[,2],0.67)))] <- "High"

# mod$type <- paste0(mod$KRAS_type,"_",mod$BMS_type)


#==========================================================================================================
# another RAS pathway 
# mutation group have higher pathway score ,However is not significant 
#===========================================================================================================
genes <- c("CAMK2B","DAB2IP","GRB2","HRAS","KRAS","LGALS1","LGALS3","NF1",
    "NRAS","PLCE1","PRKCA","PRKCB","PRKCE","PRKCZ","RABGEF1","RASA1","RASA2",
    "RASA4","RASAL1","RASGRF1","RASGRF2","RASGRP1","RASGRP2","RASGRP3",
    "RASGRP4","RIN1","RRAS","SOS1","SOS2","SYNGAP1")

source('~/software/ssGSEA/ssgseaMOD.r')
mod.generate(genes,'PID_RAS_activate',out='/public/workspace/lily/MOD_file/PID_RAS_activate.mod')
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData")
tmp <- mod.analyze2(as.matrix(dat),c("PID_RAS_activate"),'/public/workspace/lily/MOD_file/',permN=0)
tmp <- data.frame(tmp)

maffile <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_LUAD.maf',sep='\t',header=T)
KRAS_mut_sample <- unique(as.vector(maffile$Tumor_Sample_Barcode[which(maffile$Hugo_Symbol=="KRAS")]))

tmp$KRAS_type <- "non"
tmp$KRAS_type[which(rownames(tmp)%in%gsub("-",".",KRAS_mut_sample))] <- "mut"

load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData')




#=========================================================================================================
# GSE72094 data verify 
# prepare 
#=========================================================================================================
# awk -F "\t" '/^!Series_sample_id|kras_status/ {print $0}' ./GSE72094_series_matrix.txt  > KRAS_info.txt
# sed  's/kras_status: //g' ./KRAS_info.txt | sed 's/"//g' > f.txt
#=========================================================================================================
dat <- read.delim2("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_series_matrix.txt",sep="\t",comment.char="!",header=T)
info <- read.delim2("/public/workspace/lily/metastasis/data/verify/GSE72094/GPL15048_f.txt",sep="\t",header=T,comment.char="#")
# info.tmp <- info[-which(info$GeneSymbol==""),c("ID","GeneSymbol")]
# write.table(info.tmp,file="/public/workspace/lily/metastasis/data/verify/GSE72094/GPL15048_f.txt",sep="\t",col.names=T,row.names=F,quote=F)
#========================================================================
tmp <- merge(dat,info,by.x="ID_REF",by.y="ID")
tmp$ID_REF <- NULL
tmp.f <- data.frame(t(apply(tmp[,-443],1,function(x){as.numeric(as.vector(x))})))
colnames(tmp.f) <- colnames(tmp)[-443]
tmp.f$GeneSymbol <- tmp$GeneSymbol
tmp.res <- aggregate(.~GeneSymbol,data=tmp.f,FUN=median)
rownames(tmp.res) <- tmp.res$GeneSymbol
tmp.res$GeneSymbol <- NULL
saveRDS(as.matrix(tmp.res),file="/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")

#========================================================================
# check result 
# 2020-12-19
#========================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
rs <- as.data.frame(t(dat[c("MYBL2","CEBPB"),]))
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE72094/KRAS_info.txt",header=T))
ann <- ann[-1,,drop=F]
colnames(ann) <- "KRAS_type"
res.final <- merge(rs,ann,by="row.names")

#=======================
# plot result 
# density is not ok 
# use CEBPB/MYBL2 to calculate 
library(ggplot2)
library(reshape2)
#ggplot(res.final, aes(x=MYBL2, fill=KRAS_type)) + geom_density()+theme_bw()
res.final$relative_Exp <- res.final$CEBPB/res.final$MYBL2
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE72094_Fig3_mut.pdf",useDingbats=F)
boxplot(relative_Exp~KRAS_type,data=res.final,FUN=median,outline=F,ylim=c(1,2.8))
legend("topright",legend=paste0("P=",wilcox.test(relative_Exp~KRAS_type,data=res.final)$p.value))
dev.off()


#========================================================================
# KRAS up and down pathway score 
#========================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- as.data.frame(mod)
mod.f <- cbind(mod,as.data.frame(t(dat[c("MYBL2","CEBPB"),])))
mod.f$activate <- mod.f$HALLMARK_KRAS_SIGNALING_UP_norm -mod.f$HALLMARK_KRAS_SIGNALING_DN_norm
mod.f$down <- mod.f$HALLMARK_KRAS_SIGNALING_DN_norm -mod.f$HALLMARK_KRAS_SIGNALING_UP_norm
mod.f$type.act <- "high"
mod.f$type.act[which(mod.f$activate<quantile(mod.f$activate,0.5))] <- "low"

mod$CCLE_ID <- rownames(mod)



#============================================================================
# use things to calculate TF activity 
# 2020-12-23
# run in R-4.0.2
#============================================================================
library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
data(dorothea_hs,package="dorothea")
#subset DoRothEA to the confidence levels A and B to include only the high quality regulons
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
regulons = dorothea_hs[which(dorothea_hs$confidence%in%c("A","B","C","D")),]
tf_activities <- run_viper(dat, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))

rs <- as.data.frame(t(tf_activities[c("MYBL2","CEBPB"),]))
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE72094/KRAS_info.txt",header=T))
ann <- ann[-1,,drop=F]
colnames(ann) <- "KRAS_type"
res.final <- merge(rs,ann,by="row.names")


#============================================================================================================
# another GSE to verfiy 
# GSE34914
#============================================================================================================
# paste ./*.txt  > cbind.txt
dat <- read.table("~/metastasis/data/verify/GSE34914/cbind.txt",header=T,sep="\t")
dat.f <- dat[,c("Gene",colnames(dat)[grep("count",colnames(dat))])]
colnames(dat.f) <- gsub("_geneTable_start.count","",colnames(dat.f))


#=================== 
# calculate TPM 
info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",sep="\t",header=T)
colnames(info) <- c("gene","chr","start","end")
res <- merge(info,dat.f,by.x="gene",by.y="Gene")
res$length <- res$end -res$start
# to TPM 
tmp <- t(apply(res[,5:21],1,function(x){x/x[length(x)]}))*10^3
tmp.res <- apply(tmp,2,function(x){x/sum(x)})*10^6
tmp.res <- data.frame(tmp.res)
tmp.res$gene_name <- res$gene
rownames(tmp.res) <- tmp.res$gene_name
tmp.res$gene_name <- NULL
tmp.res$length <- NULL
saveRDS(tmp.res,file="~/metastasis/data/verify/GSE34914/GSE34914_expr.RDS")
#==========================================================================================================
# WT : WT <- c("LU242","LU528","LU53","LU256","LU273","LU499","LU439","LU242A")
dat <- readRDS("~/metastasis/data/verify/GSE34914/GSE34914_expr.RDS")
# calculate KRAS pathway activity
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- as.data.frame(mod)
mod.f <- cbind(mod,as.data.frame(t(dat[c("MYBL2","CEBPB"),])))
mod.f$activate <- mod.f$HALLMARK_KRAS_SIGNALING_UP_norm -mod.f$HALLMARK_KRAS_SIGNALING_DN_norm
mod.f$down <- mod.f$HALLMARK_KRAS_SIGNALING_DN_norm -mod.f$HALLMARK_KRAS_SIGNALING_UP_norm

tmp <- data.frame(t(dat[c("CEBPB","MYBL2","KRAS"),]))
tmp$KRAS_type <- "Mut"
tmp$KRAS_type[which(rownames(tmp)%in%c("LU242","LU528","LU53","LU256","LU273","LU499","LU439","LU242A"))] <- "WT"

# plot 
tmp$relative_Exp <- tmp$CEBPB/tmp$MYBL2
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE34914_Fig3_mut.pdf",useDingbats=F)
boxplot(relative_Exp~KRAS_type,data=tmp,FUN=median,outline=F)
legend("topright",legend=paste0("P=",wilcox.test(relative_Exp~KRAS_type,data=tmp)$p.value))
dev.off()





#===========================================================================================================
# other data 
# GSE63882
#
#===========================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_series_matrix.txt",sep="\t",header=T,comment.char="!")
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GPL6884.txt",sep="\t",header=T,quote="")
tmp.res <- merge(ann,dat,by.x="ID",by.y="ID_REF")
tmp.res$ID <- NULL
res.f <- aggregate(.~ILMN_Gene,data=tmp.res,FUN=median)
rownames(res.f) <- res.f$ILMN_Gene
res.f$ILMN_Gene <- NULL
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_expr.RDS")


dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_expr.RDS")
# source('~/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
# mod <- as.data.frame(mod)
# mod.f <- cbind(mod,as.data.frame(t(dat[c("MYBL2","CEBPB"),])))
# mod.f$activate <- mod.f$HALLMARK_KRAS_SIGNALING_UP_norm -mod.f$HALLMARK_KRAS_SIGNALING_DN_norm
# mod.f$down <- mod.f$HALLMARK_KRAS_SIGNALING_DN_norm -mod.f$HALLMARK_KRAS_SIGNALING_UP_norm
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_ann.txt",sep="\t",header=T)
res <- cbind(tmp,as.data.frame(t(dat[c("KRAS","CEBPB","MYBL2"),])))
colnames(res)[4] <- "KRAS_type"
#res.f <- res[which(res$Somker=="Y"),]

res.f <- res[grep("Adenocarcinoma",res$Cellline),]
res.f$relative_Exp <- res.f$CEBPB/res.f$MYBL2
# plot 
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE63882_Fig3_mut.pdf",useDingbats=F)
boxplot(relative_Exp~KRAS_type,data=res.f,FUN=median,outline=F)
legend("topright",legend=paste0("P=",wilcox.test(relative_Exp~KRAS_type,data=res.f)$p.value))
dev.off()






#============================================================================================================
# KRAS kncokdown data
# GSE15326
#============================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE15326/GSE15326_series_matrix.txt",sep="\t",header=T,comment.char="!")
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/GPL1261.txt",sep="\t",header=T)
tmp.res <- merge(tmp,dat,by.x="ID",by.y="ID_REF")

tmp.res$ENTREZ_GENE_ID <- NULL
tmp.res$ID <- NULL
res.f <- aggregate(.~Gene.Symbol,data=tmp.res,FUN=median)

res.final <- res.f[-1,]
rownames(res.final) <- res.final$Gene.Symbol
res.final$Gene.Symbol <- NULL
tmp.dat  <- apply(res.final,2,function(x){as.numeric(as.vector(x))})
rownames(tmp.dat) <- rownames(res.final)
tmp.dat.f <- tmp.dat[-grep("///",rownames(tmp.dat)),]
saveRDS(tmp.dat.f,file="/public/workspace/lily/metastasis/data/verify/GSE15326/GSE15326_expr.RDS")

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE15326/GSE15326_expr.RDS")
tmp <- data.frame(t(dat[c("Kras","Cebpb","Mybl2"),]))
tmp$relative_Exp <- tmp$Cebpb/tmp$Mybl2
tmp$type <- "Control"
tmp$type[c(1,2,9,10)] <- "shKras"
tmp$type[c(3,4)] <- "shWt1"
tmp.f <- tmp[which(tmp$type%in%c("Control","shKras")),]

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE15326_Fig3_mut.pdf",useDingbats=F)
boxplot(relative_Exp~type,data=tmp.f,FUN=median,outline=F)
legend("topright",legend=paste0("P=",wilcox.test(relative_Exp~type,data=tmp.f)$p.value))
dev.off()




apply(dat[c("Kras","Cebpb","Mybl2"),],1,function(x){c(median(x[c(1,2,9,10)]),median(x[c(5:8)]))})
















#===============================================================================================================
# GDSC analysis 
# 2020-12-24
#===============================================================================================================
RAF_pathway <- c("ARAF","BRAF","BRAP","CALM1","CAMK2A","CAMK2B",
    "CAMK2D","CAMK2G","HRAS","JAK2","KRAS","KSR1","MAP2K1","MAP2K2",
    "MAP3K11","MARK3","MRAS","NRAS","PHB","PPP1CB","PPP1CC","PPP2CA",
    "PPP2CB","PPP2R1A","PPP2R1B","PPP2R5A","PPP2R5B","PPP2R5C",
    "PPP2R5D","PPP2R5E","RAF1","SHOC2","SRC","YWHAB")

source('~/software/ssGSEA/ssgseaMOD.r')
mod.generate(RAF_pathway,'RAF_activate',out='/public/workspace/lily/MOD_file/RAF_activate.mod')

# TCGA verify 
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData")
tmp <-mod.analyze2(as.matrix(dat),c("RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
tmp <- data.frame(tmp)
maffile <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_LUAD.maf',sep='\t',header=T)
KRAS_mut_sample <- unique(as.vector(maffile$Tumor_Sample_Barcode[which(maffile$Hugo_Symbol=="KRAS")]))
tmp$KRAS_type <- "non"
tmp$KRAS_type[which(rownames(tmp)%in%gsub("-",".",KRAS_mut_sample))] <- "mut"
load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData')
tmp$BMS <- luad_mod[,2]
cor.test(tmp$RAF_activate_norm,tmp$BMS)

#===========================================================================================================
# Cell line Verify
# 81 cell line CCLE
sampleinfo <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/8mutation.txt",sep="\t",header=T)
expr <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
# filter cell line 
expr.f <- expr[,which(colnames(expr)%in%sampleinfo$CCLE_ID)]

source('~/software/ssGSEA/ssgseaMOD.r')
tmp <- mod.analyze2(as.matrix(expr.f),c("RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
tmp <- as.data.frame(tmp)
dat <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/8mut_cox_info.txt",sep="\t",header=T)
tmp$CCLE_ID <- rownames(tmp)
tmp$BMS <- mod[,2]


#===========================================================================================================
# 400 + LUAD samples 
# GSE72094
#===========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
tmp <- mod.analyze2(as.matrix(dat),c("RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
tmp.res <- cbind(tmp,mod)
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE72094/KRAS_info.txt",header=T))
ann <- ann[-1,,drop=F]
colnames(ann) <- "KRAS_type"
tmp.res <- cbind(data.frame(tmp.res),ann)
cor.test(tmp.res$RAF_activate_norm,tmp.res$BMS_test_norm)





#==========================================================================================================
# GSE63882
# 40+ cell line
#==========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
tmp <- mod.analyze2(as.matrix(dat),c("RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_ann.txt",sep="\t",header=T)
colnames(ann)[4] <- "KRAS_type"
tmp.res <- cbind(tmp,mod,ann)
tmp.res.f <- tmp.res[grep("Adenocarcinoma",tmp.res$Cellline),]
#res.f <- res[which(res$Somker=="Y"),]

res.f <- res[grep("Adenocarcinoma",res$Cellline),]
res.f$relative_Exp <- res.f$CEBPB/res.f$MYBL2



#============================================================================================================
# GDSC data 
# 2020-12-24
#============================================================================================================
library(readxl)

# GDSC1 
dat <- data.frame(read_xlsx("/public/workspace/zhangtt/GDSC1_fitted_dose_response_25Feb20.xlsx"))
luad.dat <- dat[which(dat$TCGA_DESC=="LUAD"),]
saveRDS(luad.dat,file="/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC1_LUAD_ann.RDS")

# GDSC2 
dat <- data.frame(read_xlsx("/public/workspace/zhangtt/GDSC2_fitted_dose_response_25Feb20.xlsx"))
luad.dat <- dat[which(dat$TCGA_DESC=="LUAD"),]
saveRDS(luad.dat,file="/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC2_LUAD_ann.RDS")


# cell line details
dat <- data.frame(read_xlsx("/public/workspace/zhangtt/Cell_Lines_Details.xlsx"))
luad.dat <- dat[which(dat$Cancer.Type...matching.TCGA.label=="LUAD"),] 
saveRDS(luad.dat,file="/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_ann.RDS")


# cell line expression 
dat <- read.table("/public/workspace/lily/Lung2Brain/inte7/GDSC/Cell_line_RMA_proc_basalExp.txt",sep="\t",header=T)
dat$GENE_title <- NULL
dat.f <- aggregate(.~GENE_SYMBOLS,data=dat,FUN=median)
dat.f <- dat.f[-1,]
rownames(dat.f) <- dat.f$GENE_SYMBOLS
dat.f$GENE_SYMBOLS <- NULL
saveRDS(dat.f,file="/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_expr.RDS")






#=================================================================================================
# run BMS 
#=============================
cellline.ann <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_ann.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_expr.RDS")
dat.f <- dat[,which(colnames(dat)%in%paste0("DATA.",cellline.ann$COSMIC.identifier))] #  just use LUAD cell line expression data 

# calculate BMS score 
#=============================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)

#=============================
# GDSC1 data 
gdsc1 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC1_LUAD_ann.RDS")
mod$ID <- gsub("DATA\\.","",rownames(mod))
res <- merge(gdsc1,mod,by.x="COSMIC_ID",by.y="ID")

res.f <- matrix(ncol=3)
drug <- unique(res$DRUG_NAME)
for(i in 1:length(drug)){
    tmp <- res[which(res$DRUG_NAME==drug[i]),]
    tmp.res <- c(drug[i],cor.test(tmp$BMS_test_norm,tmp$LN_IC50)$estimate,cor.test(tmp$BMS_test_norm,tmp$LN_IC50)$p.value)
    res.f <- rbind(res.f,tmp.res)
}
res.f <- res.f[-1,]
res.f <- data.frame(res.f)
rownames(res.f) <- res.f[,1]
res.f$V1 <- NULL
res.final <- apply(res.f,2,function(x){as.numeric(as.vector(x))})
rownames(res.final) <- rownames(res.f)
colnames(res.final) <- c("correlation","pvalue")
res.final <- data.frame(res.final)
res.final$log2p <- -log2(res.final$pvalue)
# rs <- res.final[which(res.final$pvalue<0.05),]
res.final$text <- ""
res.final$text[which(res.final$log2p>5)] <- rownames(res.final)[which(res.final$log2p>5)]


#=========================================================
# ggplot 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GDSC1_IC50_drug.pdf",useDingbats=F)
ggplot(res.final,aes(x=correlation,y=log2p,color=log2p))+ geom_point() +  geom_text(aes(label = text), size = 3) +
    scale_colour_gradientn(colours=c("#007cc0","#ffb310","#ed1c24")) +
    geom_hline(yintercept=5 ,linetype=4) +
    theme_classic()
dev.off()



#============================
# GDSC2 
#============================

cellline.ann <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_ann.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_expr.RDS")
dat.f <- dat[,which(colnames(dat)%in%paste0("DATA.",cellline.ann$COSMIC.identifier))] #  just use LUAD cell line expression data 

# calculate BMS score 
#=============================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
gdsc2 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC2_LUAD_ann.RDS")
mod$ID <- gsub("DATA\\.","",rownames(mod))
res <- merge(gdsc2,mod,by.x="COSMIC_ID",by.y="ID")

res.f <- matrix(ncol=3)
drug <- unique(res$DRUG_NAME)
for(i in 1:length(drug)){
    tmp <- res[which(res$DRUG_NAME==drug[i]),]
    tmp.res <- c(drug[i],cor.test(tmp$BMS_test_norm,tmp$LN_IC50)$estimate,cor.test(tmp$BMS_test_norm,tmp$LN_IC50)$p.value)
    res.f <- rbind(res.f,tmp.res)
}
res.f <- res.f[-1,]
res.f <- data.frame(res.f)
rownames(res.f) <- res.f[,1]
res.f$V1 <- NULL
res.final <- apply(res.f,2,function(x){as.numeric(as.vector(x))})
rownames(res.final) <- rownames(res.f)
colnames(res.final) <- c("correlation","pvalue")
res.final <- data.frame(res.final)
res.final$log2p <- -log2(res.final$pvalue)
# rs <- res.final[which(res.final$pvalue<0.05),]
res.final$text <- ""
res.final$text[which(res.final$log2p>5)] <- rownames(res.final)[which(res.final$log2p>5)]



pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GDSC2_IC50_drug.pdf",useDingbats=F)
ggplot(res.final,aes(x=correlation,y=log2p,color=log2p))+ geom_point() +  geom_text(aes(label = text), size = 3) +
    scale_colour_gradientn(colours=c("#007cc0","#ffb310","#ed1c24")) +
    geom_hline(yintercept=5 ,linetype=4) +
    theme_classic()
dev.off()






#==========================================================================================================================
# MEK inhbitor Verifiy #
# GSE110397
#==========================================================================================================================
tmp.dat <- read.table("~/metastasis/data/verify/GSE110397/GSE110397_NSCLC_FPKM.txt",sep="\t",header=T)
dat <- tmp.dat[,c("ensembl_gene_id","A01_S97_L002_R1_001","A02_S98_L002_R1_001","A05_S101_L002_R1_001","A06_S102_L002_R1_001",
    "A09_S105_L002_R1_001","A10_S106_L002_R1_001","A13_S109_L002_R1_001","A14_S110_L002_R1_001")]
ann <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",sep="\t",header=T)
tmp.res <- merge(dat,ann[,c(1,2),drop=F],by.x="ensembl_gene_id",by.y="Gene.stable.ID")
tmp.res$ensembl_gene_id <- NULL
res <- aggregate(.~Gene.name,data=tmp.res,FUN=median)
rownames(res) <- res$Gene.name
res$Gene.name <- NULL
saveRDS(res,file="/public/workspace/lily/metastasis/data/verify/GSE110397/GSE110397_expr.RDS")

#========================================================
# calculate BMS score 
#========================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE110397/GSE110397_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
mod$type <- c("Control","Drug","Control","Drug","Control","Drug","Control","Drug")
mod$sample <- c(rep("A549",4),rep("H2030",4))

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE110397_MEK_inhibitor.pdf",useDingbats=F)
boxplot(BMS_test_norm~type,data=mod[which(mod$sample=="A549"),],FUN=median)
boxplot(BMS_test_norm~type,data=mod[which(mod$sample=="H2030"),],FUN=median)
dev.off()











#========================================================================================================
# another GSE data 
# GSE 79235
#=======================================================
ann <- read.table("~/metastasis/data/verify/GPL17077.txt",sep="\t",header=T)
tmp.dat <- read.table("~/metastasis/data/verify/GSE79235/GSE79235_series_matrix.txt",sep="\t",header=T)
tmp.res <- merge(ann,tmp.dat,by.x="ID",by.y="ID_REF")
tmp.res$ID <- NULL
res <- aggregate(.~GENE_SYMBOL,data=tmp.res,FUN=median)
res.f <- res[-1:28,]
rownames(res.f) <- res.f$GENE_SYMBOL
res.f$GENE_SYMBOL <- NULL
saveRDS(res.f,file="~/metastasis/data/verify/GSE79235/GSE79235_expr.RDS")


dat <- readRDS("~/metastasis/data/verify/GSE79235/GSE79235_expr.RDS")
dat.f <- dat[,5:8]
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE79235_MEK_inhibitor.pdf",useDingbats=F)
boxplot(mod[1:2,2],mod[3:4,2],names=c("Control","Treat"))
dev.off()
















#==============================================================================================================
# add some result into Figure K 
# cell percentage 
# 2020-12-28
#==============================================================================================================
library(Seurat)
library(monocle)
cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
table(cds$type.TF,cds$type_group) -> tmp











# Found a KRAS KD data 
####### Wed Jan 20 15:59:59 CST 2021
#==================================================================================================================
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/KRAS_KD/counts.txt",header=T,sep="\t")
# load annotation info 
tmp.ann <- read.table("/public/workspace/lily/REF/mmGRCm38.genelength.txt",header=T,sep="\t")
tmp.ann$length <- tmp.ann$end -tmp.ann$start
tmp.res <- merge(tmp.dat,tmp.ann[,c(5,6,7)],by.x="Gene",by.y="gene_id")

# transform to TPM 
tmp <- t(apply(tmp.res[,c(2:13,15)],1,function(x){x/x[length(x)]}))*10^3
tmp.f <- apply(tmp,2,function(x){x/sum(x)})*10^6
tmp.f <- data.frame(tmp.f)
tmp.f$gene_name <- tmp.res$gene_name
final.res <- aggregate(.~gene_name,data=tmp.f,FUN=median)
rownames(final.res) <- final.res$gene_name
final.res$gene_name <- NULL
final.res$length <- NULL


res.f <- final.res[c("Kras","Cebpb","Mybl2"),]

KO <- c("SRR5099551","SRR5099553","SRR5099556","SRR5099559","SRR5099560","SRR5099561")
NC <- c("SRR5099550","SRR5099552","SRR5099554","SRR5099555","SRR5099557","SRR5099558")

KO <- c(2,4,7,10,11,12)
NC <- c(1,3,5,6,8,9)














######################################################################################################################
# analysis new LCBM data 
# scRNA
#=====================================================================================================================
library(Seurat)
tmp <- Read10X(data.dir = "/public/workspace/lily/Mutiple_LB/D0927/D0927/outs/filtered_feature_bc_matrix")
dat<- CreateSeuratObject(counts = tmp,  project = "D0927",min.cells = 3, min.features = 200)
VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dat = subset(x=dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25)
dat = NormalizeData(object = dat)
dat <- FindVariableFeatures(object = dat)
all.genes <- rownames(x = dat)
dat <- ScaleData(object = dat, features = all.genes)
# Run PCA
dat <- RunPCA(object = dat, features = VariableFeatures(object = dat),npcs = 50)
dat <- FindNeighbors(dat, dims = 1:20)
dat <- FindClusters(dat)
# the tutorial says should use the same dims in Findcluster function
dat <- RunTSNE(object = dat,dims = 1:20)
saveRDS(dat,file="/public/workspace/lily/Lung2Brain/inte7/Data/multiple_LB/D0927.RDS")












#===================================================================================================
# 2021-2-18
# analysis if KRAS is change in current and non-current samples 
# change a KRAS signal pathway
# 	Gene Set: KRAS.LUNG_UP.V1_DN
# change another KRAS pathway activaty signature 

# mod.generate(KRAS.lung.up,'KRAS.Lung.Up',out='/public/workspace/lily/MOD_file/KRAS.Lung.Up.mod')
mod.generate(gene,'RAS.pathway.act',out='/public/workspace/lily/MOD_file/RAS.pathway.act.mod')
#===================================================================================================
# 1. GSE9971 sample 
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE9971/GSE9971_series_matrix.txt",comment.char="!",header=T)
ann <- read.delim2("~/metastasis/data/verify/GPL96-57554.txt",comment.char="#",sep="\t")

tmp.dat <- merge(dat,ann,by.x="ID_REF",by.y="ID")
tmp.dat.f <- tmp.dat[-grep("///",tmp.dat$Gene.Symbol),]
tmp.dat.f$ID_REF <- NULL
tmp.res <- aggregate(.~Gene.Symbol,data=tmp.dat.f,FUN=median)
tmp.res.f <- tmp.res[-1,]
rownames(tmp.res.f) <- tmp.res.f$Gene.Symbol
tmp.res.f$Gene.Symbol <- NULL
saveRDS(tmp.res.f,file="/public/workspace/lily/metastasis/data/verify/GSE9971/GSE9971_expr.RDS")

#==================================================================================================
tmp <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE9971/GSE9971_expr.RDS")
dat <- apply(tmp,2,function(x){as.numeric(as.vector(x))})
rownames(dat)<- rownames(tmp)
saveRDS(dat,file="/public/workspace/lily/metastasis/data/verify/GSE9971/GSE9971_expr.RDS")

# check KRAS expression value
boxplot(dat["KRAS",1:11],dat["KRAS",12:27])
wilcox.test(dat["KRAS",1:11],dat["KRAS",12:27])

# calculate RAS pathway activity
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE9971/GSE9971_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("RAS.pathway.act","BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod

####################################################################################################
# 2. GSE7880
#==================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE7880/GSE7880_series_matrix.txt",sep="\t",header=T,comment.char="!")
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GPL201-30390.txt",sep="\t",header=T)


tmp.dat <- merge(dat,ann,by.x="ID_REF",by.y="ID")
tmp.dat.f <- tmp.dat[-grep("///",tmp.dat$Gene.Symbol),]
tmp.dat.f$ID_REF <- NULL
tmp.res <- aggregate(.~Gene.Symbol,data=tmp.dat.f,FUN=median)
tmp.res.f <- tmp.res[-1,]
rownames(tmp.res.f) <- tmp.res.f$Gene.Symbol
tmp.res.f$Gene.Symbol <- NULL
saveRDS(tmp.res.f,file="/public/workspace/lily/metastasis/data/verify/GSE7880/GSE7880_expr.RDS")

tmp <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE7880/GSE7880_expr.RDS")
dat <- apply(tmp,2,function(x){as.numeric(as.vector(x))})
rownames(dat)<- rownames(tmp)
saveRDS(dat,file="/public/workspace/lily/metastasis/data/verify/GSE7880/GSE7880_expr.RDS")


#===================================================================================================
# just use ademocarcinoma
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE7880/GSE7880_expr.RDS")
dat.f <- dat[,1:25]
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- data.frame(mod)

mod$activity <- mod$HALLMARK_KRAS_SIGNALING_UP_norm - mod$HALLMARK_KRAS_SIGNALING_DN_norm
wilcox.test((mod[1:10,"activity"]),(mod[11:25,"activity"]))
boxplot((mod[1:10,"activity"]),(mod[11:25,"activity"]))










#####################################################################################################
# 3. GSE8894
#====================================================================================================
# sed -n '/Series_sample_id/p' ./GSE8894_series_matrix.txt >samplename.txt
# sed -n '/Sample_characteristics_ch1/p' ./GSE8894_series_matrix.txt > sampleinfo.txt
#====================================================================================================
tmp <- read.delim2("/public/workspace/lily/metastasis/data/verify/GSE8894/sampleinfo.txt",sep="\t",header=F)
name <- read.delim2("/public/workspace/lily/metastasis/data/verify/GSE8894/samplename.txt",sep="\t",header=F)
tmp <- tmp[,-1]
colnames(tmp) <- strsplit(as.vector(name$V2)," ")[[1]]

options(stringsAsFactors=F)
ann <- data.frame(apply(tmp,1,function(x){
    sapply(as.list(as.vector(x)),function(y){strsplit(y,": ")[[1]][2]})
})) 

colnames(ann) <- c("Recurrent","Age","Gender","Type","RFS","Other")
rownames(ann) <- colnames(tmp)

# change some sample info 
ann["GSM225774",] <- c("0",NA,"Female","Adenocarcinoma",46.96667,NA)
ann["GSM225827",] <- c("1","60","Male","squamous cell carcinoma","1.5333",NA)
ann["GSM225779",] <- c("0","66","Male","squamous cell carcinoma","62.83333",NA)
ann["GSM225844",] <- c("1","67","Male","squamous cell carcinoma","6.5",NA)
ann["GSM225864",] <- c("1",NA,"Male","Adenocarcinoma","4.46667",NA)

saveRDS(ann,file="/public/workspace/lily/metastasis/data/verify/GSE8894/ann.RDS")

#=========================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_series_matrix.txt",sep="\t",header=T,comment.char="!")
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GPL570.txt",sep="\t",header=T)


tmp.dat <- merge(dat,ann,by.x="ID_REF",by.y="ID")
tmp.dat.f <- tmp.dat[-grep("///",tmp.dat$Gene.Symbol),]
tmp.dat.f$ID_REF <- NULL
tmp.res <- aggregate(.~Gene.Symbol,data=tmp.dat.f,FUN=median)
tmp.res.f <- tmp.res[-c(1:25),]
rownames(tmp.res.f) <- tmp.res.f$Gene.Symbol
tmp.res.f$Gene.Symbol <- NULL
saveRDS(tmp.res.f,file="/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")

tmp <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")
dat <- apply(tmp,2,function(x){as.numeric(as.vector(x))})
rownames(dat)<- rownames(tmp)
saveRDS(dat,file="/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")

# check 
#==========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/GSE8894_expr.RDS")
ann <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE8894/ann.RDS")
ann.f <- ann[which(ann$Type=="Adenocarcinoma"),]
dat.f <- dat[,which(colnames(dat)%in%rownames(ann.f))]
source('~/software/ssGSEA/ssgseaMOD.r')

mod <- mod.analyze2(as.matrix(dat.f),c("RAS.pathway.act","BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
# mod$activity <- mod$KRAS.Lung.Up_norm - mod$KRAS.Lung.Dn_norm

#mod <- mod.analyze2(as.matrix(dat.f),c('HALLMARK_KRAS_SIGNALING_DN',"HALLMARK_KRAS_SIGNALING_UP"),'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
# mod <- data.frame(mod)
# mod$activity <- mod$HALLMARK_KRAS_SIGNALING_UP_norm - mod$HALLMARK_KRAS_SIGNALING_DN_norm

mod$Recurrent <- ann.f$Recurrent

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE8894_KRAS_pathway_BMS.pdf",useDingbats=F)
boxplot(RAS.pathway.act_norm~Recurrent,data=mod,FUN=median)
boxplot(BMS_test_norm~Recurrent,data=mod,FUN=median)

library(ggExtra)
library(ggplot2)
library(ggpubr)

p<-ggplot(mod,aes(x=BMS_test_norm,y=RAS.pathway.act_norm)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("RAS pathway activity Score") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()






#====================================================================================================================
# 2021-2-20
# RAS pathway to check CEBPB and MYBL2



# 1. calculate activity and Expression correlation 
#====================================================================================================================
# Run in R -4.0.2

library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
data(dorothea_hs,package="dorothea")
#subset DoRothEA to the confidence levels A and B to include only the high quality regulons
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_expr.RDS")
regulons = dorothea_hs[which(dorothea_hs$confidence%in%c("A","B","C","D")),]
tf_activities <- run_viper(dat, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))

rs <- as.data.frame(t(tf_activities[c("MYBL2","CEBPB"),]))

# calculate RAS pathway activity
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("RAS.pathway.act"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- as.data.frame(mod)
# mod.f <- cbind(mod,as.data.frame(t(dat[c("MYBL2","CEBPB"),])))
# mod.f$activate <- mod.f$HALLMARK_KRAS_SIGNALING_UP_norm -mod.f$HALLMARK_KRAS_SIGNALING_DN_norm
# mod.f$down <- mod.f$HALLMARK_KRAS_SIGNALING_DN_norm -mod.f$HALLMARK_KRAS_SIGNALING_UP_norm
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_ann.txt",sep="\t",header=T)
res <- cbind(tmp,as.data.frame(t(dat[c("KRAS","CEBPB","MYBL2"),])))
res$RAS.score <- mod$RAS.pathway.act_norm
colnames(res)[4] <- "KRAS_type"
#res.f <- res[which(res$Somker=="Y"),]

res.f <- res[grep("Adenocarcinoma",res$Cellline),]
res.f$relative_Exp <- res.f$CEBPB/res.f$MYBL2

res$RAS.type <- "unknow"
res$RAS.type[which(res$RAS.score>quantile(res$RAS.score,0.5))] <- "high"
res$RAS.type[which(res$RAS.score<quantile(res$RAS.score,0.5))] <- "low"





####################################################################################################################
# GSE14108 28 LCBM sample Bulk RNA seq 
#===================================================================================================================
library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
data(dorothea_hs,package="dorothea")
#subset DoRothEA to the confidence levels A and B to include only the high quality regulons
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
regulons = dorothea_hs[which(dorothea_hs$confidence%in%c("A","B","C","D")),]
tf_activities <- run_viper(dat, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))
saveRDS(t(tf_activities),file="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_tf.RDS")

# calculata BMS score 
#==================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
tf_activities <- data.frame(readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_tf.RDS"))

tf_activities$BMS <- mod$BMS_test_norm
tmp.res <- t(apply(tf_activities,2,function(x){
    tmp <- cor.test(x,tf_activities$BMS)
    c(tmp$p.value,tmp$estimate)
}))












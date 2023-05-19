
#======================================================================================
# immunothearapy for BMS 
# 2020-8-27
#======================================================================================




#========================================================================================================
# MSK nature genetics
# 
#========================================================================================================

mskdat <- read.delim("~/MSK_LUAD/TMB/data_mutations_mskcc.txt",sep="\t",header=T)
mskdat$MAF <- mskdat$t_alt_count/(mskdat$t_ref_count+mskdat$t_alt_count) # caculaye maf 
mskdat.f <- mskdat[,c(colnames(mskdat)[1:20],"MAF")] #filter some information
mskdat.f$gene <- mskdat.f$Hugo_Symbol
#======================================== load gene BMS score 
#========================================================================================================
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(mskdat.f,res.vaf,by="gene")

#========================================================================================================
res.f <- aggregate(BMS_test_norm~Tumor_Sample_Barcode,data=res,FUN=median)
colnames(res.f)[1] <- "SAMPLE_ID"
# choose high and low risk samples 
dat.high <- res.f$SAMPLE_ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$SAMPLE_ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]

#===============================================================================================================================================================================clinical dat
mskclin <- read.delim("~/MSK_LUAD/TMB/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin[which(mskclin$CANCER_TYPE_DETAILED=="Lung Adenocarcinoma"),] -> luad_sample # all LUAD patients
#==================================================================================
sample.f <- luad_sample[which(luad_sample$SAMPLE_ID %in% c(as.vector(dat.high),as.vector(dat.low))),]

#====================survival dat
msksurv <- read.delim("~/MSK_LUAD/TMB/data_clinical_patient.txt",sep="\t",comment.char="#")
msksurv.f <- merge(msksurv,sample.f,by="PATIENT_ID")
msksurv.ff <- merge(msksurv.f,res.f,by="SAMPLE_ID")

msksurv.ff$type <- "BMS_low"
msksurv.ff$type[which(msksurv.ff$SAMPLE_ID%in%as.vector(dat.high))] <- "BMS_high"
msksurv.ff$OS_status <- 0 #change type info 
msksurv.ff$OS_status[which(msksurv.ff$OS_STATUS=="DECEASED")] <- 1
#====================
# run surv analysis 
library(survminer)
library(survival)

table(msksurv.ff$type)
surv <- Surv(msksurv.ff$OS_MONTHS,msksurv.ff$OS_status)
km <- survfit(surv~type,data=msksurv.ff)
ggsurvplot(km,pval=T)

#=================================================================
# add some filter paramater 
#=================================================================
# 1 msksurv.ff$DRUG_TYPE
msksurv.fff <- msksurv.ff[which(msksurv.ff$DRUG_TYPE=="PD-1/PDL-1"),]
table(msksurv.fff$type)
surv <- Surv(msksurv.fff$OS_MONTHS,msksurv.fff$OS_status)
km <- survfit(surv~type,data=msksurv.fff)
ggsurvplot(km,pval=T) 

# 2 GENE PANEL 
msksurv.fff <- msksurv.ff[which(msksurv.ff$GENE_PANEL=="IMPACT468"),]   # IMPACT341 IMPACT410 IMPACT468
table(msksurv.fff$type)
surv <- Surv(msksurv.fff$OS_MONTHS,msksurv.fff$OS_status)
km <- survfit(surv~type,data=msksurv.fff)
ggsurvplot(km,pval=T) 


# 3 SAMPLE_TYPE
msksurv.fff <- msksurv.ff[which(msksurv.ff$SAMPLE_TYPE=="Primary"),]   # IMPACT341 IMPACT410 IMPACT468
table(msksurv.fff$type)
surv <- Surv(msksurv.fff$OS_MONTHS,msksurv.fff$OS_status)
km <- survfit(surv~type,data=msksurv.fff)
ggsurvplot(km,pval=T) 


#====================================
msksurv.fff <- msksurv.ff[which(msksurv.ff$DRUG_TYPE=="PD-1/PDL-1"&msksurv.ff$GENE_PANEL=="IMPACT468"&msksurv.ff$SAMPLE_TYPE=="Metastasis"),]
table(msksurv.fff$type)
surv <- Surv(msksurv.fff$OS_MONTHS,msksurv.fff$OS_status)
km <- survfit(surv~type,data=msksurv.fff)
ggsurvplot(km,pval=T) 


#========================================================================================================
# check gene panel 
# not very difference 
#========================================================================================================
msk410 <- read.table("./MSK_LUAD/TMB/MSK410gene.txt")
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
colnames(msk410)[1] <- "gene"
msk410.f <- merge(msk410,res.vaf,by='gene')

load("/public/workspace/lily/MSK_LUAD/TMB/gene.rda")
msk468 <- as.data.frame(a)
colnames(msk468)[1] <- "gene"
msk468.f <- merge(msk468,res.vaf,by="gene")












#==========================================================================================================
#
# GSE126044
#==========================================================================================================
#==========================================================================================================
# calculate ssgsea
# 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE126044/TPM.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("LUADgene","LUSCgene"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- as.data.frame(mod)

mod2 <- mod.analyze2(as.matrix(dat),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod2 <- as.data.frame(mod2)


mod.f <- mod2[which(rownames(mod2)%in%c("Dis_01","Dis_11","Dis_15","Dis_04","Dis_16","Dis_07","Dis_10","Dis_05","Dis_08")),]
mod.f$type <- "non-responder"
mod.f$type[which(rownames(mod.f)%in%c("Dis_02","Dis_15","Dis_04","Dis_17","Dis_10"))] <- "responder"
aggregate(BMS_test_norm~type,data=mod.f,FUN=median)



#==========================================================================================================
# filter LUSC 
# when filter LUSC ,responder BMS norm is higher
# 
#==========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE126044/TPM.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat.f <- dat[,which(colnames(dat)%in%c("Dis_01","Dis_11","Dis_15","Dis_04","Dis_16","Dis_07","Dis_10","Dis_05","Dis_08"))]
mod.f <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod.f <- as.data.frame(mod.f)
mod.f$type <- "non-responder"
mod.f$type[which(rownames(mod.f)%in%c("Dis_02","Dis_15","Dis_04","Dis_17","Dis_10"))] <- "responder"
aggregate(BMS_test_norm~type,data=mod.f,FUN=median)








































#===================================================================================================================================================================================
# MSK 240 NSCLC 2018 IMPACT
#===================================================================================================================================================================================

# ways 1 : do not use gene ,but use gene BMS score to caculate sample BMS score 
#=============================================================================================
# msk mutation data
mskdat <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_mutations_mskcc.txt",sep="\t",header=T)
mskdat$MAF <- mskdat$t_alt_count/(mskdat$t_ref_count+mskdat$t_alt_count) # caculaye maf 
mskdat.f <- mskdat[,c(colnames(mskdat)[1:20],"MAF")] #filter some information
mskdat.f$Hugo_Symbol -> mskdat.f$gene
#=============================================================================================
# gene BMS score 
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(mskdat.f,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Tumor_Sample_Barcode,data=res,FUN=mean)
dat.high <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]))
dat.low <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]))


#================================================================================================
# msk clinical dat
mskclin <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin.f <- mskclin[which(mskclin$ONCOTREE_CODE=="LUAD"&mskclin$GENE_PANEL=="IMPACT410"),] # IMPACT341  IMPACT410  IMPACT468   &mskclin$GENE_PANEL=="IMPACT341"
mskclin.f$type <- "unsure"
mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.high)] <- "High"
mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.low)] <- "Low"
if(length(which(mskclin.f$type=="unsure"))>0){
	mskclin.f <- mskclin.f[-which(mskclin.f$type=="unsure"),]
}

#====================survival dat
msksurv <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_patient.txt",sep="\t",comment.char="#")
msksurv.f <- merge(msksurv,mskclin.f,by="PATIENT_ID")



#===========================
# check TMB 
mskclin <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin.f <- mskclin[which(mskclin$ONCOTREE_CODE=="LUAD"),] 
mskclin.f$type <- "unsure"
mskclin.f$type[which(mskclin.f$MUTATION_RATE>quantile(mskclin.f$MUTATION_RATE,0.5))] <- "High"
mskclin.f$type[which(mskclin.f$MUTATION_RATE<quantile(mskclin.f$MUTATION_RATE,0.5))] <- "Low"
msksurv <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_patient.txt",sep="\t",comment.char="#")
msksurv.f <- merge(msksurv,mskclin.f,by="PATIENT_ID")




#============================================================================================================================
# TMB high with BMS 
# 
#============================================================================================================================
mskclin <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin.f <- mskclin[which(mskclin$ONCOTREE_CODE=="LUAD"),] 
mskclin.f <- mskclin.f[which(mskclin.f$MUTATION_RATE>quantile(mskclin.f$MUTATION_RATE,0.5)),]

mskdat <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_mutations_mskcc.txt",sep="\t",header=T)
mskdat$MAF <- mskdat$t_alt_count/(mskdat$t_ref_count+mskdat$t_alt_count) # caculaye maf 
mskdat.f <- mskdat[,c(colnames(mskdat)[1:20],"MAF")] #filter some information
mskdat.f$Hugo_Symbol -> mskdat.f$gene
mskdat.f <- mskdat.f[which(mskdat.f$Tumor_Sample_Barcode%in%mskclin.f$SAMPLE_ID),]

res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(mskdat.f,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Tumor_Sample_Barcode,data=res,FUN=mean)
dat.high <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]))
dat.low <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]))


mskclin.f <- mskclin[which(mskclin$ONCOTREE_CODE=="LUAD"&mskclin$GENE_PANEL=="IMPACT410"),] # IMPACT341  IMPACT410  IMPACT468   &mskclin$GENE_PANEL=="IMPACT341"
mskclin.f$type <- "unsure"
mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.high)] <- "High"
mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.low)] <- "Low"
if(length(which(mskclin.f$type=="unsure"))>0){
	mskclin.f <- mskclin.f[-which(mskclin.f$type=="unsure"),]
}

#====================survival dat
msksurv <- read.delim("/public/workspace/lily/MSK_LUAD/MSK240NSCLC/data_clinical_patient.txt",sep="\t",comment.char="#")
msksurv.f <- merge(msksurv,mskclin.f,by="PATIENT_ID")












#====================
# run surv analysis 
library(survminer)
library(survival)

msksurv.f$PFS <- as.numeric(as.vector(sapply(strsplit(as.vector(msksurv.f$PFS_STATUS),":"),function(x){x[1]})))
table(msksurv.f$type)
surv <- Surv(msksurv.f$PFS_MONTHS,msksurv.f$PFS)
km <- survfit(surv~type,data=msksurv.f)
ggsurvplot(km,pval=T)


msksurv.ff <- msksurv.f[which(msksurv.f$TREATMENT_TYPE=="Monotherapy"),]
surv <- Surv(msksurv.ff$PFS_MONTHS,msksurv.ff$PFS)
km <- survfit(surv~type,data=msksurv.ff)

table(msksurv.f$type,msksurv.f$DURABLE_CLINICAL_BENEFIT)
table(msksurv.ff$type,msksurv.ff$DURABLE_CLINICAL_BENEFIT)

surv <- Surv(msksurv.ff$PFS_MONTHS,msksurv.ff$PFS)
km <- survfit(surv~type,data=msksurv.ff)
ggsurvplot(km,pval=T)








#========================================================================================================================================
# MSK 970 NSCLC 2017
# the data do not have many immunotherapy data
# so do not use 
#========================================================================================================================================
# ways 1 : do not use gene ,but use gene BMS score to caculate sample BMS score 
#=============================================================================================
# # msk mutation data
# mskdat <- read.delim("/public/workspace/lily/MSK_LUAD/MSK970NSCLC/data_mutations_mskcc.txt",sep="\t",header=T,comment.char="#")
# mskdat$MAF <- mskdat$t_alt_count/(mskdat$t_ref_count+mskdat$t_alt_count) # caculaye maf 
# mskdat.f <- mskdat[,c(colnames(mskdat)[1:20],"MAF")] #filter some information
# mskdat.f$Hugo_Symbol -> mskdat.f$gene
# #=============================================================================================
# # gene BMS score 
# res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
# res <- merge(mskdat.f,res.vaf,by="gene")
# res.f <- aggregate(BMS_test_norm~Tumor_Sample_Barcode,data=res,FUN=mean)
# dat.high <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]))
# dat.low <- as.vector(unique(res.f$Tumor_Sample_Barcode[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]))
# #================================================================================================
# # msk clinical dat
# mskclin <- read.delim("/public/workspace/lily/MSK_LUAD/MSK970NSCLC/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
# mskclin.f <- mskclin[which(mskclin$ONCOTREE_CODE=="LUAD"),] # IMPACT341  IMPACT410  IMPACT468   &mskclin$GENE_PANEL=="IMPACT341"
# mskclin.f$type <- "unsure"
# mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.high)] <- "High"
# mskclin.f$type[which(mskclin.f$SAMPLE_ID%in%dat.low)] <- "Low"
# if(length(which(mskclin.f$type=="unsure"))>0){
# 	mskclin.f <- mskclin.f[-which(mskclin.f$type=="unsure"),]
# }

# #====================survival dat
# msksurv <- read.delim("/public/workspace/lily/MSK_LUAD/MSK970NSCLC/data_clinical_patient.txt",sep="\t",comment.char="#")
# msksurv.f <- merge(msksurv,mskclin.f,by="PATIENT_ID")
# msksurv.ff <- msksurv.f[which(msksurv.f$IMMUNE_TREATMENT=="YES"&msksurv.f$SAMPLE_TYPE=="Primary"),]  #  Primary  Metastasis

# table(msksurv.ff$type,msksurv.ff$DURABLE_CLINICAL_BENEFIT)



#=======================================================================================================================================
# GSE135222
# NSCLC sample : maybe have LUSC samples 
#=======================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/60sample_mutation.txt",sep="\t",header=T)
sapply(strsplit(gsub("\\(",":",as.vector(dat$Gene)),":"),function(x){x[1]}) -> dat$gene
luadsample <- c("NSCLC378","NSCLC1203","NSCLC1401","NSCLC825","NSCLC1528","NSCLC1017","NSCLC1358","NSCLC1510","NSCLC1164","NSCLC1145","NSCLC1619","NSCLC1809","NSCLC1873")
TMB <- table(dat$Sample.ID)



#===================================================
# surv info 
#===================================================
sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.txt",sep="\t",header=T)
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(dat,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Sample.ID,data=res,FUN=mean)
res.f$sample <- paste0("NSCLC",res.f$Sample.ID)
res.f <- res.f[which(res.f$sample%in%luadsample),]
dat.high <- res.f$Sample.ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$Sample.ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]
#====================================================
# add info  ID is number ,so can use this 
sample.info[which(sample.info$Patient.ID%in%c(dat.high,dat.low)),] -> surv.info
surv.info$type <- "low"
surv.info$type[which(surv.info$Patient.ID%in%c(dat.high))] <- "high"

#=============================

library(survminer)
library(survival)

surv <- Surv(surv.info$PFS,surv.info$PFS_event)
km <- survfit(surv~type,data=surv.info)



#==========================
# check TMB 
#
#==========================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/60sample_mutation.txt",sep="\t",header=T)
sapply(strsplit(gsub("\\(",":",as.vector(dat$Gene)),":"),function(x){x[1]}) -> dat$gene
luadsample <- c("NSCLC378","NSCLC1203","NSCLC1401","NSCLC825","NSCLC1528","NSCLC1017","NSCLC1358","NSCLC1510","NSCLC1164","NSCLC1145","NSCLC1619","NSCLC1809","NSCLC1873")
TMB <- data.frame(table(dat$Sample.ID))
colnames(TMB)[1] <- "Patient.ID"

sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.txt",sep="\t",header=T)
res.f <- merge(TMB,sample.info,by="Patient.ID")
res.f$sample <- paste0("NSCLC",res.f$Patient.ID)
res.f <- res.f[which(res.f$sample%in%luadsample),]
dat.high <- res.f$Patient.ID[which(res.f$Freq>quantile(res.f$Freq,0.5))]
dat.low <- res.f$Patient.ID[which(res.f$Freq<quantile(res.f$Freq,0.5))]

sample.info[which(sample.info$Patient.ID%in%c(as.numeric(as.vector(dat.high)),as.numeric(as.vector(dat.low)))),] -> surv.info
surv.info$type <- "low"
surv.info$type[which(surv.info$Patient.ID%in%as.numeric(as.vector(dat.high)))] <- "high"

#=============================

library(survminer)
library(survival)

surv <- Surv(surv.info$PFS,surv.info$PFS_event)
km <- survfit(surv~type,data=surv.info)









#===================================================================================================================================
# cancer cell 
#===================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/mutation.txt",sep="\t",header=T)
sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/sample.txt",header=T,sep="\t")
sample.info.f <- sample.info[which(sample.info$Histology=="non-squamous"),] # non-squamous 
#================== gene score 
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
dat$gene <- dat$gene_name
res <- merge(dat,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Patient_ID,data=res,FUN=mean)
dat.high <- res.f$Patient_ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$Patient_ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]
# filter 
dat.high.f <- dat.high
dat.low.f <- dat.low
# ============================
info.surv <- sample.info.f[which(sample.info.f$Patient.ID%in%c(dat.high.f,dat.low.f)),c("Patient.ID","PFS..months.","PFS...0.censor..1.event.","Best.Overall.Response","Clinical.benefit..DCB...durable.clinical.benefit..NDB...no.durable.benefit")]
info.surv$type <- "low"
info.surv$type[which(info.surv$Patient.ID%in%dat.high.f)] <- "high"
colnames(info.surv) <- c("Patient.ID","PFS.time","PFS","Response","DCB","type")
#=============================
library(survminer)
library(survival)
surv <- Surv(info.surv$PFS.time,info.surv$PFS)
km <- survfit(surv~type,data=info.surv)
ggsurvplot(km,pval=T)




#==================================Check TMB 
#================================================================
TMB.high <- sample.info.f$Patient.ID[which(sample.info.f$Nonsynonymous.tumor.mutation.burden>quantile(sample.info.f$Nonsynonymous.tumor.mutation.burden,0.5))]
TMB.low <- sample.info.f$Patient.ID[which(sample.info.f$Nonsynonymous.tumor.mutation.burden<quantile(sample.info.f$Nonsynonymous.tumor.mutation.burden,0.5))]


dat.high.f <- TMB.high
dat.low.f <- TMB.low
# ============================
info.surv <- sample.info.f[which(sample.info.f$Patient.ID%in%c(dat.high.f,dat.low.f)),c("Patient.ID","PFS..months.","PFS...0.censor..1.event.","Best.Overall.Response","Clinical.benefit..DCB...durable.clinical.benefit..NDB...no.durable.benefit")]
info.surv$type <- "low"
info.surv$type[which(info.surv$Patient.ID%in%dat.high.f)] <- "high"
colnames(info.surv) <- c("Patient.ID","PFS.time","PFS","Response","DCB","type")
#=============================
library(survminer)
library(survival)
surv <- Surv(info.surv$PFS.time,info.surv$PFS)
km <- survfit(surv~type,data=info.surv)
ggsurvplot(km,pval=T)




#================================TMB 
#================================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/mutation.txt",sep="\t",header=T)
sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/sample.txt",header=T,sep="\t")
sample.info.f <- sample.info[which(sample.info$Histology=="non-squamous"),]
sample.info.f <- sample.info.f[which(sample.info.f$Nonsynonymous.tumor.mutation.burden<quantile(sample.info.f$Nonsynonymous.tumor.mutation.burden,0.5)),]
dat.f <- dat[which(dat$Patient_ID%in%sample.info.f$Patient.ID),]
#====================================
# BMS score 
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
dat.f$gene <- dat.f$gene_name
res <- merge(dat.f,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Patient_ID,data=res,FUN=mean)
dat.high <- res.f$Patient_ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$Patient_ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]

# filter 
dat.high.f <- dat.high
dat.low.f <- dat.low
# ============================
info.surv <- sample.info.f[which(sample.info.f$Patient.ID%in%c(dat.high.f,dat.low.f)),c("Patient.ID","PFS..months.","PFS...0.censor..1.event.","Best.Overall.Response","Clinical.benefit..DCB...durable.clinical.benefit..NDB...no.durable.benefit")]
info.surv$type <- "low"
info.surv$type[which(info.surv$Patient.ID%in%dat.high.f)] <- "high"
colnames(info.surv) <- c("Patient.ID","PFS.time","PFS","Response","DCB","type")
#=============================
library(survminer)
library(survival)
surv <- Surv(info.surv$PFS.time,info.surv$PFS)
km <- survfit(surv~type,data=info.surv)
ggsurvplot(km,pval=T)




























#================================================================================================================================================
# use LUAD and LUSC to classifiy 
# geneset to classifiy 
# CCAR genelist 
#================================================================================================================================================
LUADgene <- c("SPINK1","SFTA2","NKX2-1","CAPN8","CLDN3","TESC","LGSN","HNF1B","GOLT1A","HPN","FMO5","ACSL5","CDH15","ALDH3B1","RORC","DPP4","PRR15L","STK32A","KCNK5","ABCC6","SMPDL3B","FAM83B")
LUSCgene <- c("DST","COL7A1","SOX15","PNCK","VSNL1","DSC3","KRT16","SERPINB5","DAPL1","KRT17","PKP1","DSG3","KRT13","CLCA2","KRT6C","TP63","KRT6B","CALML3","KRT6A","KRT5")

genelist <- rep(1,length=(length(LUADgene)+length(LUSCgene)))
names(genelist) <- c(LUADgene,LUSCgene)
genelist[which(names(genelist)%in%LUSCgene)] <- (-1)

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(genelist,"LUAD_LUSC_sig",out="/public/workspace/lily/MOD_file/LUAD_LUSC_sig.mod")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(LUADgene,"LUADgene",out="/public/workspace/lily/MOD_file/LUADgene.mod")
mod.generate(LUSCgene,"LUSCgene",out="/public/workspace/lily/MOD_file/LUSCgene.mod")


dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/expr.RDS")
mod2 <- mod.analyze2(as.matrix(dat),c("LUADgene","LUSCgene"),"/public/workspace/lily/MOD_file/",permN=1000)
mod2 <- as.data.frame(mod2)

#===============================================================================================
# load data and calculate 
#===============================================================================================
luadsample <- c("NSCLC378","NSCLC1203","NSCLC1401","NSCLC825","NSCLC1528","NSCLC1017","NSCLC1358","NSCLC1510","NSCLC1164","NSCLC1145","NSCLC1619","NSCLC1809","NSCLC1873")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/expr.RDS")
dat.f <- dat[,which(colnames(dat)%in%luadsample)]
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- as.data.frame(mod)
mod.f <- mod
surv <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.RDS")
res <- merge(mod.f,surv,by="row.names")

res$type <- "median"
res$type[which(res$BMS_test_norm<quantile(res$BMS_test_norm,0.5))] <- "low"
res$type[which(res$BMS_test_norm>=quantile(res$BMS_test_norm,0.5))] <- "high"

library(survminer)
library(survival)
surv <- Surv(res$PFS,res$PFS_event)
km <- survfit(surv~type,data=res)
ggsurvplot(km,pval=T)






# # more strict value 
# #==============================================================================
# luadsample <- c("NSCLC378","NSCLC1203","NSCLC1401","NSCLC825","NSCLC1528","NSCLC1510","NSCLC1164","NSCLC1809","NSCLC1873")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/expr.RDS")
# dat.f <- dat[,which(colnames(dat)%in%luadsample)]
# mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
# mod <- as.data.frame(mod)

# mod[which(rownames(mod)%in%luadsample),] -> mod.f
# surv <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.RDS")

# res <- merge(mod.f,surv,by="row.names")

# res$type <- "median"
# res$type[which(res$BMS_test_norm<quantile(res$BMS_test_norm,0.5))] <- "low"
# res$type[which(res$BMS_test_norm>=quantile(res$BMS_test_norm,0.5))] <- "high"

# library(survminer)
# library(survival)
# surv <- Surv(res$PFS,res$PFS_event)
# km <- survfit(surv~type,data=res)
# ggsurvplot(km,pval=T)






























































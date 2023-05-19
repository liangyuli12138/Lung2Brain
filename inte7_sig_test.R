
#========================================================================================================================================================================
# llt write in 2020-8-5
# used for gene signature test analysis 
# run in pang .8
#========================================================================================================================================================================
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

#========================================================================================================================================================================
# get mod file path and modname
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
modpath <- "/public/workspace/lily/Lung2Brain/inte7/" # maybe Lung2Bain use a same file folder 
modname <- "BMS_test"
#mod.generate(modgene,modname,out=paste0(modpath,modname,".mod")) # generate mod file 
dir.create(paste0("/public/workspace/lily/Lung2Brain/Opt_sig/",modname))
resfile <- paste0("/public/workspace/lily/Lung2Brain/Opt_sig/",modname,"/")
setwd(resfile)


#====================================================================
# run ss-gsea for 3 vs 3 sample 
#====================================================================
library(Seurat)
load('/public/workspace/lily/metastasis/data/verify/GSE126548/GSE126548_rpkm.RData')
sig_mod <- mod.analyze2(as.matrix(rpkm),modname,modpath,permN=0)
save(sig_mod,file=paste0(resfile,'GSE126548_sig_mod','.rda'))

#====================================================================
# plot boxplot
#====================================================================
pdf(paste0(resfile,'sig_boxplot33','.pdf')) #boxplot for GSE126548 (3 vs. 3)
# sig1
boxplot(sig_mod[c(1:3),2],sig_mod[c(4:6),2],names=c('BM-','BM+'),main=paste0(modname))
legend('topright',legend=paste('p =',round(wilcox.test(sig_mod[c(1:3),2],sig_mod[c(4:6),2])$p.value,5)),bty='n')
dev.off()


#====================================================================
# run survival analysize -- run ssgsea
#====================================================================
# load data 
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/TCGA_LUSC_data.RData') # lusc_dat
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # dat 
dat -> luad_dat
rm(dat) #change variable names
# run ssgsea
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
luad_mod <- mod.analyze2(as.matrix(luad_dat),modname,modpath,permN=0)
lusc_mod <- mod.analyze2(as.matrix(lusc_dat),modname,modpath,permN=0)

save(luad_mod,file=paste0(resfile,'LUAD_sig_mod','.rda'))
save(lusc_mod,file=paste0(resfile,'LUSC_sig_mod','.rda'))

#====================================================================
# run survival analysize -- plot 
#====================================================================
#==================================
# survival
#==================================
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
ann -> luad_ann
rm(ann)
#
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/TCGA_LUSC_ann.RData')
ann -> lusc_ann
rm(ann)
#
# load(paste0(resfile,'LUAD_sig_mod','.rda'))
# load(paste0(resfile,'LUSC_sig_mod','.rda'))

############ sig1 #################
# add type info 
# LUAD five quantile
luad_ann$type_five <- rep("M",nrow(luad_ann))
luad_ann[which(luad_mod[,2]>unname(quantile(luad_mod[,2],0.8))),"type_five"] <- "H"
luad_ann[which(luad_mod[,2]<unname(quantile(luad_mod[,2],0.2))),"type_five"] <- "L"
# LUAD half
luad_ann$type_two <- rep("L",nrow(luad_ann))
luad_ann[which(luad_mod[,2]>unname(quantile(luad_mod[,2],0.5))),"type_two"] <- "H"
# LUAD three quantile
luad_ann$type_three <- rep("M",nrow(luad_ann))
luad_ann[which(luad_mod[,2]>unname(quantile(luad_mod[,2],0.67))),"type_three"] <- "H"
luad_ann[which(luad_mod[,2]<unname(quantile(luad_mod[,2],0.33))),"type_three"] <- "L"

# LUSC five quantile 
lusc_ann$type_five <- rep("M",nrow(lusc_ann))
lusc_ann[which(lusc_mod[,2]>unname(quantile(lusc_mod[,2],0.8))),"type_five"] <- "H"
lusc_ann[which(lusc_mod[,2]<unname(quantile(lusc_mod[,2],0.2))),"type_five"] <- "L"
# LUSC half
lusc_ann$type_two <- rep("L",nrow(lusc_ann))
lusc_ann[which(lusc_mod[,2]>unname(quantile(lusc_mod[,2],0.5))),"type_two"] <- "H"
# LUSC three quantile 
lusc_ann$type_three <- rep("M",nrow(lusc_ann))
lusc_ann[which(lusc_mod[,2]>unname(quantile(lusc_mod[,2],0.67))),"type_three"] <- "H"
lusc_ann[which(lusc_mod[,2]<unname(quantile(lusc_mod[,2],0.33))),"type_three"] <- "L"

############ plot #################
pdf(paste0("./",'sig_surv_LUAD','.pdf'),useDingbats=F)
#==================================
# plot OS time
# sig LUAD type : half 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survfit(surv~type_two,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_half",xlab=" Overall survival (months)")+labs(title="LUAD_sig_half")
# sig LUAD type : five 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survfit(surv~type_five,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_five",xlab=" Overall survival (months)")+labs(title="LUAD_sig_five")
# sig LUAD type : three 
surv <- Surv(luad_ann$OS.time,luad_ann$OS)
km <- survfit(surv~type_three,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" Overall survival (months)")+labs(title="LUAD_sig_three")

# plot RFS time
# sig1 LUAD type : half 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survfit(surv~type_two,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_half",xlab=" RF survival (months)")+labs(title="LUAD_sig_half")
# sig1 LUAD type1 : five 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survfit(surv~type_five,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_five",xlab=" RF survival (months)")+labs(title="LUAD_sig_five")
# sig1 LUAD type2 : three 
surv <- Surv(luad_ann$RFS_time,luad_ann$RFS)
km <- survfit(surv~type_three,data=luad_ann)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',xlim=c(0,3000),xscale='d_m',break.time.by =304.37,
	title="LUAD_sig_three",xlab=" RF survival (months)")+labs(title="LUAD_sig_three")
dev.off()

#######################################################################################################################################################################################
pdf(paste0(resfile,'sig_surv_LUSC','.pdf'))
#==================================
# plot OS time
# sig1 LUSC type : half 
surv <- Surv(lusc_ann$OS.time,lusc_ann$OS)
km <- survfit(surv~type_two,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_half",xlab=" Overall survival (months)")+labs(title="LUSC_sig_half")
# sig1 LUSC type : five 
surv <- Surv(lusc_ann$OS.time,lusc_ann$OS)
km <- survfit(surv~type_five,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_five",xlab=" Overall survival (months)")+labs(title="LUSC_sig_five")
# sig1 LUSC type : three 
surv <- Surv(lusc_ann$OS.time,lusc_ann$OS)
km <- survfit(surv~type_three,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_three",xlab=" Overall survival (months)")+labs(title="LUSC_sig_three")

# plot RFS time
# sig1 LUSC type : half 
surv <- Surv(lusc_ann$RFS.time,lusc_ann$RFS)
km <- survfit(surv~type_two,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_half",xlab=" RF survival (months)")+labs(title="LUSC_sig_half")
# sig1 LUSC type : five 
surv <- Surv(lusc_ann$RFS.time,lusc_ann$RFS)
km <- survfit(surv~type_five,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_five",xlab=" RF survival (months)")+labs(title="LUSC_sig_five")
# sig1 LUSC type : three 
surv <- Surv(lusc_ann$RFS.time,lusc_ann$RFS)
km <- survfit(surv~type_three,data=lusc_ann)
ggsurvplot(km,pval=T,palette=c('#F8766D','#00BFC4','grey'),surv.median.line='hv',xlim=c(0,2000),xscale='d_m',break.time.by =304.37,
	title="LUSC_sig_three",xlab=" RF survival (months)")+labs(title="LUSC_sig_three")
dev.off()





#====================================================================
# TCGA metasis analysize
#====================================================================
library(ggplot2)
library(ggpubr)


pdf(paste0(resfile,'sig_TCGA_Metastasis','.pdf'))
##############sig1#################
# LUAD
luad_clin <- read.table('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt',sep='\t',header=T)
luad_clin <- luad_clin[,c(1,102:104)]
my_comparisons <- list(c("M1", "M0"))
# add type info 
clin <- data.frame(sample=gsub('-','.',luad_clin$sampleID),type_M=luad_clin$pathologic_M,type_N=luad_clin$pathologic_N)
mod <- data.frame(sample=rownames(luad_mod),score=luad_mod[,2])
res <- merge(clin,mod,by='sample')
res$type <- as.vector(res$type_M)
res$type[intersect(which(res$type_M=='M0'),which(res$type_N=="N1"))] <- 'N1'
res$type[which(res$type_M=='M1')] <- 'M1'
res[which(res$type%in%c('M1',"M0","N1")),] -> res_f
# plot
ggplot(data=res_f,aes(x=type,y=score,fill=type))+
    stat_boxplot(geom = "errorbar",width=0.25)+
	geom_boxplot(outlier.shape = NA,outlier.size=1)+
	stat_compare_means(comparisons=my_comparisons,method.args = list(alternative = "greater"))+
	theme_classic()+labs(title='sig_LUAD_greater')

# LUSC
lusc_clin <- read.table('/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/LUSC_clinicalMatrix',sep='\t',header=T)
lusc_clin <- lusc_clin[,c(1,74:76)]
my_comparisons <- list(c("M1", "M0"))
# add type info 
clin <- data.frame(sample=gsub('-','.',lusc_clin$sampleID),type_M=lusc_clin$pathologic_M,type_N=lusc_clin$pathologic_N)
mod <- data.frame(sample=rownames(lusc_mod),score=lusc_mod[,2])
res <- merge(clin,mod,by='sample')
res$type <- as.vector(res$type_M)
res$type[intersect(which(res$type_M=='M0'),which(res$type_N=="N1"))] <- 'N1'
res$type[which(res$type_M=='M1')] <- 'M1'
res[which(res$type%in%c('M1',"M0","N1")),] -> res_f
# plot
ggplot(data=res_f,aes(x=type,y=score,fill=type))+
    stat_boxplot(geom = "errorbar",width=0.25)+
	geom_boxplot(outlier.shape = NA,outlier.size=1)+
	stat_compare_means(comparisons=my_comparisons,method.args = list(alternative = "greater"))+
	theme_classic()+labs(title='sig_LUSC_greater')
dev.off()




#===========================================================
# GSE131907 validate
#===========================================================
library(Seurat)
library(ggplot2)
library(ggpubr)

tumor <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE131907/LungT_advance_noPTPRC.rds")
mod_lungae <- mod.analyze2(as.matrix(tumor[['RNA']]@data),modname,modpath,permN=0)
save(mod_lungae,file=paste0(resfile,"GSE131907_sig_mod.rda"))

mod_lungae <- data.frame(mod_lungae)
mod_lungae$type <- tumor$orig.ident
pdf(paste0(resfile,"GSE131907_sig_boxplot.pdf"))
my_comparisons =list(c("LungTumor_advanced","LungTumor_early"))
ggplot(data = mod_lungae,aes(x=type,y=mod_lungae[,2],fill=type))+
	geom_boxplot(outlier.shape = NA,outlier.size=1)+
	#stat_boxplot(geom = "errorbar",width=0.25)+
	stat_compare_means(comparisons=my_comparisons)+ # Add pairwise 
	labs(x="Type",y="BMS",title=modname)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))
dev.off()




#====================================================================================
# TCGA test
#====================================================================================

# load data for intraInvasion and localInvasion
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/TCGA_LUSC_data.RData') # lusc_dat
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # dat 
dat -> luad_dat
rm(dat) #change variable names
# run ssgsea
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
luad_vasion_mod <- mod.analyze2(as.matrix(luad_dat),c("localInvasion_cp",'intraInvasion_cp'),"/public/workspace/lily/MOD_file/",permN=0)
lusc_vasion_mod <- mod.analyze2(as.matrix(lusc_dat),c("localInvasion_cp",'intraInvasion_cp'),"/public/workspace/lily/MOD_file/",permN=0)

#============
# BMS score
load(paste0(resfile,'LUSC_sig_mod','.rda'))
load(paste0(resfile,'LUAD_sig_mod','.rda'))


res_luad <- merge(luad_mod,luad_vasion_mod,by='row.names')
res_lusc <- merge(lusc_mod,lusc_vasion_mod,by="row.names")

#================================================================
# result plot 
#================================================================

pdf(paste0(resfile,"TCGA_vasion.pdf"))
ggplot(data =res_luad,aes(x =res_luad[,3],y =res_luad[,7])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="local",title='TCGA_LUAD_local')
  
ggplot(data =res_luad,aes(x =res_luad[,3],y =res_luad[,8])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="extra",title='TCGA_LUAD_extra')

ggplot(data =res_lusc,aes(x =res_lusc[,3],y =res_lusc[,7])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="local",title='TCGA_LUSC_local')

ggplot(data =res_lusc,aes(x =res_lusc[,3],y =res_lusc[,8])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="extra",title='TCGA_LUSC_extra')

 
dev.off()



# load data for blood and cell migration 
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/TCGA_LUSC_data.RData') # lusc_dat
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # dat 
dat -> luad_dat
rm(dat) #change variable names
# run ssgsea
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
luad_function_mod <- mod.analyze2(as.matrix(luad_dat),c("blood_signature",'cell_migration_sig'),"/public/workspace/lily/MOD_file/",permN=0)
lusc_function_mod <- mod.analyze2(as.matrix(lusc_dat),c("blood_signature",'cell_migration_sig'),"/public/workspace/lily/MOD_file/",permN=0)

#============
# BMS score
load(paste0(resfile,'LUSC_sig_mod','.rda'))
load(paste0(resfile,'LUAD_sig_mod','.rda'))


res_luad <- merge(luad_mod,luad_function_mod,by='row.names')
res_lusc <- merge(lusc_mod,lusc_function_mod,by="row.names")


#========================================================================================================================================================================================================
# plot result 
#========================================================================================================================================================================================================
# plot BMS functions 

pdf(paste0(resfile,"TCGA_function_BMS.pdf"))
ggplot(data =res_luad,aes(x =res_luad[,3],y =res_luad[,7])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="local",title='TCGA_LUAD_angiogenesis')
  
ggplot(data =res_luad,aes(x =res_luad[,3],y =res_luad[,8])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="extra",title='TCGA_LUAD_cell_migration')

ggplot(data =res_lusc,aes(x =res_lusc[,3],y =res_lusc[,7])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="local",title='TCGA_LUSC_angiogenesis')

ggplot(data =res_lusc,aes(x =res_lusc[,3],y =res_lusc[,8])) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="extra",title='TCGA_LUSC_cell_migration')

 
dev.off()






#=============================================================
# KRAS mutation 
#=============================================================
# KRAS maftools plot
library(ggplot2)
library(maftools)

load(paste0(resfile,'LUAD_sig_mod','.rda'))
luad_mod <- luad_mod[order(luad_mod[,2],decreasing=T),]
a <- luad_mod[,2]
nM <- luad_mod[which(luad_mod[,2]<unname(quantile(a,0.33))),]
M <- luad_mod[which(luad_mod[,2]>unname(quantile(a,0.67))),]
rownames(nM) <- gsub("\\.",'-',rownames(nM))
rownames(M) <- gsub("\\.","-",rownames(M))

maffile <- read.table('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/mutation_LUAD.maf',sep='\t',header=T)
maffile[which(maffile$Tumor_Sample_Barcode%in%rownames(nM)),]-> maf_nM
maffile[which(maffile$Tumor_Sample_Barcode%in%rownames(M)),]-> maf_M
maf_n <-  read.maf(maf_nM)
maf_y <- read.maf(maf_M)

pdf(paste0(resfile,"LUAD_KRAS_mut.pdf"))
lollipopPlot(maf_n,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
printCount=T,showMutationRate = TRUE,showDomainLabel=F)

lollipopPlot(maf_y,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
printCount=T,showMutationRate = TRUE,showDomainLabel=F)
dev.off()




#=================================================================
# stemness id 
modpath <- "/public/workspace/lily/Lung2Brain/inte6/" # maybe Lung2Bain use a same file folder 
modname <- "gene_test2"
resfile <- paste0("/public/workspace/lily/Lung2Brain/Opt_sig/",modname,"/")

load(paste0(resfile,'LUAD_sig_mod','.rda'))
load(paste0(resfile,'LUSC_sig_mod','.rda'))

#=================================================================
load("/public/workspace/lily/metastasis/data/verify/STEM_id/LUAD_mRNA.RData")
rownames(luad_mrna) <- luad_mrna$TCGAlong.id
tmp <- merge(luad_mod,luad_mrna,by='row.names')
cor.test(tmp$gene_test7_norm,tmp$mRNAsi)



















































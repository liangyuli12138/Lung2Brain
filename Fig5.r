
# this program is used to calculate Figure 5 which calculate Mutation and BMS . therapy
# 2021-5-10
#================================================================================================
# 0. claulate BMS score associate mutation 
# LUAD <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
# source('~/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(LUAD),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
# mod <- data.frame(mod)
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")


#===============================================================================================
luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
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
# filter some mutation and rank by varation 
res.f <- luad_res1[which(luad_res1$p_value<0.01),]
res.ff <- res.f[which(res.f$percent>0.15),]
res.ff$varation <- (res.ff$mut_mean-res.ff$no_mut_mean)
res.ff <- res.ff[order(res.ff$varation,decreasing=T),]
res.ff$gene <- factor(res.ff$gene,levels=res.ff$gene)



#==========================================================================================================================
# 
library(ggplot2)
# ggplot(res.ff,aes(x=varation,y=percent)) + geom_point()
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/Mutation_KRAS.pdf",useDingbats=F)
ggplot(res.ff,aes(x=gene,y=varation,fill= percent)) + geom_bar(stat="identity",position="dodge") + theme_bw()
dev.off()





# 0.1 maftools check mutation point 
#==========================================================================================================================
#########KRAS##########
#library(g3viz)
library(ggplot2)
library(maftools)
#########LUAD-C1##########
luad_mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
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


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/KRAS_maftools.pdf",useDingbats=F)
lollipopPlot(maf_n,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
printCount=T,showMutationRate = TRUE,showDomainLabel=F)

lollipopPlot(maf_y,gene = "KRAS",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,printCount=T,
showMutationRate = TRUE,showDomainLabel=F) 
dev.off()




# pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TP53_maftools.pdf",useDingbats=F)
# lollipopPlot(maf_n,gene = "TP53",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,
# printCount=T,showMutationRate = TRUE,showDomainLabel=F)

# lollipopPlot(maf_y,gene = "TP53",labelPos = "all",AACol = "Amino_Acid_Change",labPosSize = 0.9,printCount=T,repel=T,
# showMutationRate = TRUE,showDomainLabel=F) 
# dev.off()








# 1. use CCLE data to verfiy [KRAS Mutation sample ] 
#=======================================================================================================================
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_KRAS_mut_sample.txt",sep="\t",header=F)
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_ann.txt",header=T,sep="\t")
tmp.f <- tmp.dat[which(tmp.dat$tcga_code=="LUAD"&tmp.dat$Pathology=="primary"),]
# tmp.f <- tmp.dat[which(tmp.dat$tcga_code=="LUAD"),]
tmp.f$KRAS_Mut <- "N"
tmp.f$KRAS_Mut[which(tmp.f$depMapID%in%tmp$V1)] <- "Y"


dat <- readRDS("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RPKM.RDS")
dat.f <- dat[,which(colnames(dat)%in%tmp.f$CCLE_ID)]

# calculate BMS score 
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
rownames(tmp.f) <- tmp.f$CCLE_ID
tmp.f <- tmp.f[rownames(mod),]
mod$KRAS_Mut <- tmp.f$KRAS_Mut
boxplot(BMS_update_norm~KRAS_Mut,data=mod,FUN=median)

# have trend but not significant



# 2. use cell line data to verify 
# have trend but not significant 
#==========================================================================================================
# GSE63882
# 40+ cell line
#==========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_expr.RDS")
ann <- read.table("/public/workspace/lily/metastasis/data/verify/GSE63882/GSE63882_ann.txt",sep="\t",header=T)
colnames(ann)[4] <- "KRAS_type"
rownames(ann) <- ann$GSM
dat <- dat[,rownames(ann)]
dat.f <- dat[,grep("Adenocarcinoma",ann$Cellline)]
ann.f <- ann[grep("Adenocarcinoma",ann$Cellline),]
#res.f <- res[which(res$Somker=="Y"),]
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update","RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
mod$KRAS_type <- ann.f$KRAS_type




# 3. use BULK data to verfiy ,signifcant
#===========================================================================================================
# 400 + LUAD samples 
# GSE72094
#===========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
source('~/software/ssGSEA/ssgseaMOD.r')
# tmp <- mod.analyze2(as.matrix(dat),c("RAF_activate"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- mod.analyze2(as.matrix(dat),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE72094/KRAS_info.txt",header=T))
ann <- ann[-1,,drop=F]
colnames(ann) <- "KRAS_type"
all(rownames(mod)==rownames(ann))
tmp.res <- cbind(mod,ann)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE72094_verfiy_KRAS.pdf",useDingbats=F)
boxplot(BMS_update_norm~KRAS_type,data=tmp.res,FUN=median,outline=F)
dev.off()



# 4.GSE15326
# trend but not significant 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE15326/GSE15326_expr.RDS")
dat.f <- dat[,c(1,2,9,10,5,6,7,8)]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update.m"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
wilcox.test(mod[1:4,2],mod[5:8,2])



# 5. GDSC 
#=============================
cellline.ann <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_ann.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_expr.RDS")
dat.f <- dat[,which(colnames(dat)%in%paste0("DATA.",cellline.ann$COSMIC.identifier))] #  just use LUAD cell line expression data 

# calculate BMS score 
#=============================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
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
    tmp.res <- c(drug[i],cor.test(tmp$BMS_update_norm,tmp$LN_IC50)$estimate,cor.test(tmp$BMS_update_norm,tmp$LN_IC50)$p.value)
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


# plot result 
# 
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GDSC1.pdf",useDingbats=F)
ggplot(res.final,aes(x=correlation,y=log2p,color=log2p))+ geom_point(size=3) +  geom_text(aes(label = text), size = 3) +
    scale_colour_gradientn(colours=c("#007cc0","#ffb310","#ed1c24")) +
    geom_hline(yintercept=5 ,linetype=4) +
    theme_classic()
dev.off()






########## GDSC2 
# 2021-5-19 
# GDSC2 result is more good
cellline.ann <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_ann.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC_Cellline_expr.RDS")
dat.f <- dat[,which(colnames(dat)%in%paste0("DATA.",cellline.ann$COSMIC.identifier))] #  just use LUAD cell line expression data 

# calculate BMS score 
#=============================
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),'/public/workspace/lily/MOD_file/',permN=0)
mod <- data.frame(mod)
gdsc2 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC2_LUAD_ann.RDS")
mod$ID <- gsub("DATA\\.","",rownames(mod))
res <- merge(gdsc2,mod,by.x="COSMIC_ID",by.y="ID")

res.f <- matrix(ncol=3)
drug <- unique(res$DRUG_NAME)
for(i in 1:length(drug)){
    tmp <- res[which(res$DRUG_NAME==drug[i]),]
    tmp.res <- c(drug[i],cor.test(tmp$BMS_update_norm,tmp$LN_IC50)$estimate,cor.test(tmp$BMS_update_norm,tmp$LN_IC50)$p.value)
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

############################################################################################
# plot result 
# 
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GDSC2.pdf",useDingbats=F)
ggplot(res.final,aes(x=correlation,y=log2p,color=log2p))+ geom_point(size=3) +  geom_text(aes(label = text), size = 3) +
    scale_colour_gradientn(colours=c("#007cc0","#ffb310","#ed1c24")) +
    geom_hline(yintercept=5 ,linetype=4) +
    theme_classic()
dev.off()








# 6. MEK inhibitor result 
# 2021-5-19
#===========================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE79235/GSE79235_expr.RDS")
dat.f <- dat[,5:8]
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE79235_MEK_inhibitor.pdf",useDingbats=F)
barplot(c(median(mod[1:2,2]),median(mod[3:4,2])),names=c("Control","Treat"),ylim=c(0,1))
dev.off()













































































###########################################################################################################################################
# 2021-7-1
# use this program to analysis difference of three subtype tumor cells in LCBM 
#==========================================================================================================================================
# 0. pheatmap show DEGs of three group 
# dose percentage have differences in LCBMs? 
# 1. metabolism difference in three group 
# 2. three subtype in DTC to macrometastasis 
# 3. differnece for immunotherapy 
# 4. 





# 1. pheatmap show DEGs 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$group)
gene <- FindAllMarkers(dat,assay="RNA",logfc.threshold=0,only.pos=T)
gene <- gene[order(gene$avg_logFC,decreasing=T),]
gene.f <- gene[which(gene$p_val_adj<0.05),]
gene.f <- gene.f[-grep("^RPL|^RPS",gene.f$gene),]
# get genes 
gene.f <- gene.f[order(gene.f$avg_logFC,decreasing=T),]
gene.use <- c()
top <- 10
for(i in unique(gene.f$cluster)){
	tmp.gene <- gene.f[which(gene.f$cluster==i),"gene"]
	gene.use <- c(gene.use,tmp.gene[1:top])
}


# plot data @data[heatmap.gene,]
# data <- dat[["RNA"]]@data[gene.use,]
# library(pheatmap)
# ann <- dat@meta.data[,"group",drop=F]
# ann <- ann[order(ann$group),,drop=F]
# pheatmap(data[,rownames(ann)],annotation_col=ann,
#     cluster_cols=F,scale="row",show_colnames=F,color=colorRampPalette(c('steelblue',"steelblue",'white',"red",'red'))(100))
# data <- data[,rownames(ann)]

# group by cluster 
data <- dat[["RNA"]]@data[gene.use,]
data <- t(as.matrix(data))
data <- data.frame(data)
data$group <- dat$group
tmp <- aggregate(.~group,data=data,FUN=mean)
rownames(tmp) <- tmp$group
tmp$group <- NULL
tmp.res <- t(tmp)


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_sub_tumor_heatmap_DEG.pdf",useDingbats=F)
pheatmap::pheatmap(tmp.res,cluster_cols=F,scale="row",
	show_colnames=T,color=colorRampPalette(c('steelblue',"steelblue",'white',"red",'red'))(100))
dev.off()






# 2. try radioresistance signature 
# this siganture come from MSGDB 
# HP_INCREASED_SENSITIVITY_TO_IONIZING_RADIATION
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("Radio_sens","Radio_sen_test","Radio_sen_BMC","Radio_sen_CCR"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$group <- dat$group
mod$group <- factor(mod$group,levels=c("group16","group09","group47"))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/Radio_sensetivity.pdf",useDingbats=F)
boxplot(Radio_sens_norm~group,data=mod,outline=F,ylim=c(0,1))
dev.off()

# 2021-10-14 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/Radio_sensetivity_new.pdf",useDingbats=F)
boxplot(Radio_sen_test_norm~group,data=mod,outline=F,ylim=c(0,1))
dev.off()


# calculate in Bulk 
# 2021-10-14
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","Radio_sen_test"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


# plot result 
gene <- c("AR","PRKCA","RELA","SUMO1","CDK1","HDAC1","IRF1")
dat <- dat.BM[gene,]
tmp.ann <- t(sapply(gene,function(x){
    c(cor.test(mod$Brain_gene_norm,as.numeric(dat[x,]),method="spearman")$estimate,
    cor.test(mod$Lung_gene_norm,as.numeric(dat[x,]),method="spearman")$estimate,
    cor.test(mod$BMS_update_norm,as.numeric(dat[x,]),method="spearman")$estimate
    )
}))
colnames(tmp.ann) <- c("Brain","Lung","Stem")
ann_row <- data.frame(apply(tmp.ann,2,function(x){ifelse(x>0,"pos","neg")}))
ann_colors = list(
    Brain = c(pos = "steelblue", neg = "red"),
    Lung = c(pos = "steelblue", neg = "red"),
    Stem = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = colorRampPalette(c('white',"#44AC48"))(100),
    Brain_gene_norm = colorRampPalette(c('white',"#269DCB"))(100),
    BMS_update_norm = colorRampPalette(c('white',"#F6BB3D"))(100)
)

ann_colors = list(
    Brain = c(pos = "steelblue", neg = "red"),
    Lung = c(pos = "steelblue", neg = "red"),
    Stem = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = c(High = "red", Medium = "white",Low = "green"),
    Brain_gene_norm = c(High = "red", Medium = "white",Low = "green"),
    BMS_update_norm = c(High = "red", Medium = "white",Low = "green")
)

ann_col <- data.frame(mod[,c(5:7)])
ann_col <- data.frame(apply(ann_col,2,function(x){ifelse(x<quantile(x,0.33),"Low",ifelse(x>quantile(x,0.67),"High","Medium"))}))
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_sentive_gene_heatmap.pdf",useDingbats=F,height=10,width=10)
pheatmap(dat,scale="row",annotation_col=ann_col,annotation_row=ann_row,annotation_colors=ann_colors,
    color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),cellwidth=8,cellheight=15
)
dev.off()


# E-MTAB
# 2021-10-14
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Lung_gene","Brain_gene","BMS_update","Radio_sens","Radio_sen_test"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

# plot result 
gene <- c("AR","PRKCA","RELA","SUMO1","CDC2","HDAC1","IRF1")
dat <- dat.BM[gene,]
tmp.ann <- t(sapply(gene,function(x){
    c(cor.test(mod$Brain_gene_norm,as.numeric(dat[x,]))$estimate,
    cor.test(mod$Lung_gene_norm,as.numeric(dat[x,]))$estimate,
    cor.test(mod$BMS_update_norm,as.numeric(dat[x,]))$estimate
    )
}))
colnames(tmp.ann) <- c("Brain","Lung","Stem")
ann_row <- data.frame(apply(tmp.ann,2,function(x){ifelse(x>0,"pos","neg")}))
ann_colors = list(
    Brain = c(pos = "steelblue", neg = "red"),
    Lung = c(pos = "steelblue", neg = "red"),
    Stem = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = colorRampPalette(c('white',"#44AC48"))(100),
    Brain_gene_norm = colorRampPalette(c('white',"#269DCB"))(100),
    BMS_update_norm = colorRampPalette(c('white',"#F6BB3D"))(100)
)
ann_col <- data.frame(mod[,c(6:8)])
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/E_MTAB_sentive_gene_heatmap.pdf",useDingbats=F,height=10,width=10)
pheatmap(dat,scale="row",annotation_col=ann_col,annotation_row=ann_row,annotation_colors=ann_colors,
    color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),cellwidth=8,cellheight=15
)
dev.off()





# RRRS
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("RES","SEN","blood_signature"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$group <- dat$group
mod$group <- factor(mod$group,levels=c("group16","group09","group47"))
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_suntumor_RRRS_mod.RDS")





mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_suntumor_RRRS_mod.RDS")
mod$Radio <- mod$RES_norm - mod$SEN_norm
mod$group <- factor(mod$group,levels=c("group09","group16","group47"))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/Radio_response.pdf",useDingbats=F)
boxplot(Radio~group,data=mod,outline=F)
dev.off()














# 3. DTC 2 Macrometastases verify 
# this was used in Figure4 
# however these samples do not have brain metastasis,maybe not every good 
#====================================================================================================================================





#====================================================================================================================================
# 4. metabolism pathway analysis 
# one is metabolism pathway activation numbers 

#====================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/sub.Tumor.metabolism.RDS")
tmp.f <- data.frame(tmp.dat[,86:170])
tmp.f$group <- dat$group

tmp.res <- aggregate(.~group,data=tmp.f,FUN=median)
rownames(tmp.res) <- tmp.res$group
tmp.res$group <- NULL

tmp.res <- t(tmp.res)
res.f <- tmp.res[-which(rowSums(tmp.res)==0),]

# prepare for plot result 
rownames(res.f) <- gsub("_norm$","",rownames(res.f))
rownames(res.f) <- gsub("_",".",rownames(res.f))

# circle plot 
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subTumor_metabolism.pdf",useDingbats=F,height=15)
plot_Circ_heatmap(t(res.f))
dev.off()


# calculate pathway numbers
tmp <- unlist(apply(res.f,1,function(x){
    colnames(res.f)[which(x[]>0.5)]
}))





# Find the Brain specific metabolism pathway and do GSEA
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
dat@active.ident <- factor(dat$group)
gene <- FindMarkers(dat,assays="RNA",ident.1="group16",min.pct=0.05,logfc.threshold=0)
gene <- gene[order(gene$avg_logFC,decreasing=T),]
write.table(gene[,2,drop=F],file="/public/workspace/lily/Lung2Brain/Version5/Fig5/Group16_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")

gene.f <- rownames(gene[which(gene$p_val_adj<0.05&gene$avg_logFC>0.25),])
write.table(gene.f,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/Group16_DEG_tmp.txt",row.names=F,col.names=F,quote=F,sep="\t")


gene <- gene[order(gene$p_val,decreasing=F),]
gene$logp <- -log()
write.table(gene[,1,drop=F],file="/public/workspace/lily/Lung2Brain/Version5/Fig5/Group16_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")





# 2021-7-14
# GSE LCBM data show Brain signature correlation with Ketone metabolism
#=======================================================================================================================================
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("Synthesis_and_degradation_of_ketone_bodies","Lung_gene","Brain_gene","BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_Brain_ketone.pdf",useDingbats=F)
plot(mod$Brain_gene_norm,mod$Synthesis_and_degradation_of_ketone_bodies_norm,mian="Brain signature with VIM")
abline(lm(mod$Synthesis_and_degradation_of_ketone_bodies_norm~mod$Brain_gene_norm),col="red")
legend("topright",legend=paste0("rho=",0.51," pvalue=",0.007))
dev.off()



# self data, E0927 and GSE123902 data use to verify 
#############################################################################################################################################
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.rs <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("Synthesis_and_degradation_of_ketone_bodies"),"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)
mod.rs$group <- factor(dat$group,levels=c("group16","group47","group09"))

# plot result 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subTumor_ketone_metabolism.pdf",useDingbats=F)
boxplot(Synthesis_and_degradation_of_ketone_bodies_norm~group,data=mod.rs,FUN=median,outline=F,ylim=c(0,1))
#wilcox.test(mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group16")],mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group47")])
dev.off()


# E0927 data 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/E0927_subTumor_ketone_metabolism.pdf",useDingbats=F)
boxplot(Synthesis_and_degradation_of_ketone_bodies_norm~group,data=mod.rs,FUN=median,outline=F,ylim=c(0,1))
#wilcox.test(mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group16")],mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group47")])
dev.off()


# GSE123902 data 
mod.rs$group <- dat$tmp.group
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_subTumor_ketone_metabolism.pdf",useDingbats=F)
boxplot(Synthesis_and_degradation_of_ketone_bodies_norm~group,data=mod.rs,FUN=median,outline=F,ylim=c(0,1))
#wilcox.test(mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group16")],mod.rs$Synthesis_and_degradation_of_ketone_bodies_norm[which(mod.rs$group=="group47")])
dev.off()







##############################################################################################################################################################
# try to find group16 expression gene and do network
# 2021-7-8
# verify by other samples show Group16 is low gene expression 
#=============================================================================================================================================================
# 2021-7-9 
# do metabolism PCA 
# however , maybe should try to use raw data 
#=============================================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/GTEx/GTEx_brain_lung.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
mod.rs <- mod.analyze2(as.matrix(dat),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)

saveRDS(mod.rs,file="/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GTEx_Lung_Brain_metabolism_mod.RDS")

#=============================================================================================================================================================
# check result 
# GTEx result 
gtex <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/GTEx_Lung_Brain_metabolism_mod.RDS")[,1:85]
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/sub.Tumor.metabolism.RDS")[,1:85]













# Analysis Therapy induced lung cancer
# seem not vevry intersting
# 2021-7-16 
#==========================================================================================================================================================
library(Seurat)
dat <- readRDS("~/metastasis/data/verify/TKI_multiple_Lung/multiple_LB_tumor.RDS")
dat@meta.data <- dat@meta.data[,c(1:14,19:21,26,40,41,45:46)]
saveRDS(dat,file="/public/workspace/lily/metastasis/data/Cell_TKI_lung/Tumor.RDS")

#==========================================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/Cell_TKI_lung/Tumor.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update","HPSC_C5","Brain_gene","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)









###########################################################################################################################################################
# 2021-7-27 
# CSOmap to run LCBM sub tumor 
#==========================================================================================================================================================
# bytlib load languages/R-4.0.2
# load package
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo)


# function 1 
as_matrix <- function(mat){
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

# function 2
run_CSOmap <- function(labelData, TPM, sampling=0, seed=315) {
    if (sampling) {
        TPM_ <- TPM
        labelData_ <- labelData

        set.seed(seed)
        index <- unlist(
            apply(as.data.frame(table(labelData_$labels) * sampling/ncol(TPM_)), 1, function(x) {
                sample(which(labelData_$labels==x[1]), x[2])
            })
        )

        TPM <- TPM_[,index]
        labelData <- labelData_[index,]
    }

    # Calculate optimized 3D coordinates
    affinityMat = getAffinityMat(TPM, LR, verbose = T)

    coords_res = runExactTSNE_R(
    X = affinityMat,
    no_dims = 3,
    max_iter = 1000,
    verbose = T
    )
    coords = coords_res$Y
    rownames(coords) <- colnames(TPM)
    colnames(coords) <- c('x', 'y', 'z')

    # Visualization(by 3D density)
    require(dplyr)
    # arrange data
    coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

    join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
    cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

    density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
    cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

    # p_3Ddensity = plot3D(cellinfo_tbl, color_by = "density", title = "3D density")

    # Get significance
    signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
    contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)

    return(list(
        labelData=labelData,
        TPM=TPM,
        cellinfo_tbl=cellinfo_tbl,
        signif_results=signif_results,
        contribution_list=contribution_list
    ))

}


tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
TPM <- as.matrix(tmp.res$RNA@counts)

labelData <- data.frame(
    cells  = colnames(tmp.res),
    labels = as.vector(tmp.res$group)
)

results = run_CSOmap(labelData, TPM, sampling=0)


saveRDS(results, file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_subtumor_CSOmap.results.RDS")


# download to local and plot 
plot_ly(result$cellinfo_tbl, x = ~x, y = ~y, z = ~z, color = ~labels , text = ~labels,colors = c("#2BAF2B","#00A5DC","#FFC719"),size=0.3)





#=============================================================================================================
# caculate different cell type distance with sub tumor 
# 2021-8-4
#-============================================================================================================
# bytlib load languages/R-4.0.2
# load package
library(Seurat)
library(CSOmapR)
library(CSOmapR.demo)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/LCBM_celltype_type.RDS")
# load function1 and function2 
TPM <- as.matrix(dat$RNA@counts)

labelData <- data.frame(
    cells  = colnames(dat),
    labels = as.vector(dat$celltype)
)

results = run_CSOmap(labelData, TPM, sampling=0)


saveRDS(results, file="/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_subtumor_celltype_CSOmap.results.RDS")
# caculate distance 
tmp <- aggregate(.~labels,data=results$cellinfo_tbl[,c(2:5)],FUN=mean)
rownames(tmp) <- tmp$labels
tmp$labels <- NULL
# download to local and plot 
plot_ly(result$cellinfo_tbl, x = ~x, y = ~y, z = ~z, color = ~labels , text = ~labels,colors = c("#2BAF2B","#00A5DC","#FFC719"),size=0.3)














###########################################################################################################################################
# 2021-7-27 
# calculate LCBM sub tumor metabolism with GBM and Lung cancer 
#==========================================================================================================================================
# LCBM sub tumor has run metabolism 
# now run GBM and LUAD

library(Seurat)

# GBM tumor 
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")
tumor <- subset(tmp.dat,cells=which(tmp.dat$type=="malignant"))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
mod.rs <- mod.analyze2(as.matrix(tumor[["RNA"]]@data),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)
saveRDS(mod.rs,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GBM_tumor_metabolism_ssGSEA.RDS")


###### LUAD 
tumor <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LUAD_clust.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
filename <- gsub("\\.mod","",dir("/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism"))
mod.rs <- mod.analyze2(as.matrix(tumor[["RNA"]]@data),filename,"/public/workspace/zhumy/CRC2Liver/ref/mod/metabolism/",permN=0)
mod.rs <- as.data.frame(mod.rs)
saveRDS(mod.rs,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/LUAD_tumor_metabolism_ssGSEA.RDS")


#================================================================================================================================================
# get result 
dat1 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GBM_tumor_metabolism_ssGSEA.RDS")
res1 <- apply(dat1[,86:170],2,function(x){median(x)})

dat2 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/LUAD_tumor_metabolism_ssGSEA.RDS")
res2 <- apply(dat2[,86:170],2,function(x){median(x)})

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/sub.Tumor.metabolism.RDS")
dat <- dat[,86:170]
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
dat$group <- tmp.dat$group
aggregate(.~group,data=dat,FUN=median)-> tmp.res
rownames(tmp.res) <- tmp.res$group
tmp.res$group <- NULL
res <- data.frame(t(tmp.res))
all(rownames(res)==names(res1))
res$GBM <- unname(res1)
res$LUAD <- unname(res2)

saveRDS(res,file="/public/workspace/lily/Lung2Brain/Version5/Fig5/GBM.LCBM.LUAD.tumor.metabolism.RDS")


# cluster show not good result 
# try to use PCA, also do not show good result 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GBM.LCBM.LUAD.tumor.metabolism.RDS")
pheatmap::pheatmap(dat[-15,],scale="row")


















#================================================================================================================================================
# 2021-7-30
# NKX2.2 in normal brain devlopment 
#================================================================================================================================================
tmp.dat <- data.table::fread("/public/workspace/lily/metastasis/data/verify/E-MTAB-6814/E-MTAB-6814-query-results.tpms.tsv",sep="\t",header=T)
tmp.f <- as.data.frame(tmp.dat)[,c(2,grep("brain",colnames(tmp.dat)))]
colnames(tmp.f) <- sapply(strsplit(colnames(tmp.f)," "),function(x){paste0(x[length(x)],"_",x[1])})
res <- tmp.f[grep("NKX2-2",tmp.f$Name_Gene),-1]
# set order 
res.f <- as.numeric(res[,c(11,1:8,12,9,10,13:31)])
pdf("~/NKX2.2_brain_dev.pdf")
plot(res.f[1:20],type="o",main="forebrain")
plot(res.f[21:31],type="o",main="hindbrain")
dev.off()




























#####################################################################################################################################################
# 2021-8-14
# Lung like tumor analysis 
# use tumor inflamation score verify in GSE14108
#====================================================================================================================================================
# make a new mod 
# https://jitc.bmj.com/content/6/1/63
gene <- c("PSMB10","HLA-DQA1","HLA-DRB1","CMKLR1","HLA-E","NKG7","CD8A",
    "CCL5","CXCL9","CD27","CXCR6","IDO1","STAT1","TIGIT","LAG3","CD274","PDCD1LG2","CD276")


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"TIS",out="/public/workspace/lily/MOD_file/TIS.mod") # make a mod file 



dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("TIS","Lung_gene","Brain_gene","BMS_update","Angiogenesis","RES","SEN"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


tmp <- read.table("~/metastasis/data/verify/GSE14108/ImmuCellAI_abundance_result_GSE14108.txt",sep="\t",header=T)
all(rownames(mod)==tmp$X)
tmp.res <- data.frame(t(apply(tmp[,-1],2,function(x){
    c(cor.test(x,mod$Brain_gene_norm,method="spearman")$estimate,cor.test(x,mod$Brain_gene_norm,method="spearman")$p.value,
        cor.test(x,mod$Lung_gene_norm,method="spearman")$estimate,cor.test(x,mod$Lung_gene_norm,method="spearman")$p.value,
        cor.test(x,mod$BMS_update_norm,method="spearman")$estimate,cor.test(x,mod$BMS_update_norm,method="spearman")$p.value)
})))
colnames(tmp.res) <- c("B.cor","B.p","L.cor","L.p","BMS.c","BMS.p")




















#====================================================================================================================================================
# 2021-8-14
# cancer cell 6A 
# use LCBM
#====================================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
tcell <- c("CD2","CD3D","CD3E","CD3G")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(tcell,"Tcell",out="/public/workspace/lily/MOD_file/Tcell.mod") # make a mod file 

# load T cell specifc 
tcell <- as.vector(read.table("/public/workspace/lily/MOD_file/Tcell_specific.txt")[,1])
dat.tcell <- t(dat[which(rownames(dat)%in% tcell ),])
# calculate score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","Brain_gene","BMS_update","Tcell"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
res <- cbind(mod[,c(5,6,7,8)],t(dat[which(rownames(dat)%in% tcell ),]))

#calculate correlation
tmp <- t(apply(res,2,function(x){
    c(
        cor.test(x,res[,1])$estimate,cor.test(x,res[,1])$p.value,
        cor.test(x,res[,2])$estimate,cor.test(x,res[,2])$p.value,
        cor.test(x,res[,3])$estimate,cor.test(x,res[,3])$p.value,
        cor.test(x,res[,4])$estimate,cor.test(x,res[,4])$p.value
    )
}))

colnames(tmp) <- c("Lung_gene.cor","Lung_gene.pval","Brain_gene.cor","Brain_gene.pval","BMS.cor","BMS.pval","Tcell.cor","Tcell.pval")
tmp <- data.frame(tmp)









#===========================================================================================================================================
# 2021-8-15 
# make lung sig and brain sig
#===========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$group)
gene <- FindAllMarkers(dat,min.pct=0.1,logfc.threshold=0.1)
gene <- gene[order(gene$avg_logFC,decreasing=T),]
# filter RP gene
gene.f <- gene[-grep("^RP",rownames(gene)),]

# get lung.gene
lung.gene <- head(gene.f[which(gene.f$cluster=="group09"),"gene"],10)
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(lung.gene,"Lung.sig",out="/public/workspace/lily/MOD_file/Lung.sig.mod") # make a mod file 

# get brain.gene
brain.gene <- head(gene.f[which(gene.f$cluster=="group16"),"gene"],10)
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(brain.gene,"Brain.sig",out="/public/workspace/lily/MOD_file/Brain.sig.mod") # make a mod file 




dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
tcell <- c("CD2","CD3D","CD3E","CD3G")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(tcell,"Tcell",out="/public/workspace/lily/MOD_file/Tcell.mod") # make a mod file 

# load T cell specifc 
tcell <- as.vector(read.table("/public/workspace/lily/MOD_file/Tcell_specific.txt")[,1])
dat.tcell <- t(dat[which(rownames(dat)%in% tcell ),])
# calculate score 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","Brain_gene","BMS_update","TIS","BMDM_marker","MG_marker","Treg","Tcell_tumor"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
res <- cbind(mod[,c(7:12)],t(dat[which(rownames(dat)%in% tcell ),]))











#==========================================================================================================================================================
# 2021-8-16 
# try to use immunosuppression signature to analysis [get OK result]
# now try to plot results
#==========================================================================================================================================================
# immunegene <- c("ADORA2A","ARHGEF5","BTLA","CD160","CD244","CD27","CD274","CD276","CD47","CD80","CEACAM1","CTLA4",
#     "GEM","HAVCR2","ICOS","IDO1","LAG3","PDCD1","TNFSF4","VISTA","VTCN1")

# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(immunegene,"immunecheck",out="/public/workspace/lily/MOD_file/immunecheck.mod") # make a mod file 

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","RES","SEN","immunecheck","Brain_gene","BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$RRS <- mod$RES_norm-mod$SEN_norm


# plot result for therapy resistance
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_RRRS_lung.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,mod$RRS,main="Lung like tumor with RRRS",xlab="Lung cancer signature",ylab="RRS signature",pch=19)
abline(lm(mod$RRS~mod$Lung_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.70," pvalue<",0.001))
dev.off()


# plot result for immunecheckpoint

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_immunecheck_lung.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,mod$immunecheck_norm,main="Lung like tumor with immunecheckpoint",xlab="Lung cancer signature",ylab="immunocheckpoint signature",pch=19)
abline(lm(mod$immunecheck_norm~mod$Lung_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.45," pvalue=",0.017))
dev.off()


# use cibersort inflammation or immuncell AI result to do correlation 
# cibersort do not have inflatrate score, so use ImmuncellAI
source("~/software/Cibersort_R.R")
tmp <- read.table("~/metastasis/data/verify/GSE14108/ImmuCellAI_abundance_result_GSE14108.txt",sep="\t",header=T)
all(rownames(tmp$X)==rownames(mod))

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_inflammation_lung.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,tmp$InfiltrationScore,main="Lung like tumor with inflammation",xlab="Lung cancer signature",ylab="inflammation signature",pch=19)
abline(lm(tmp$InfiltrationScore~mod$Lung_gene_norm),col="red")
# cor.tes()
legend("topright",legend=paste0("rho=",0.49," pvalue=",0.01))
dev.off()





#===================================================================================================================================
# check immune checkpoint in bulk data 
# 2021-10-14
####################################################################################################################################
library(pheatmap)
BM.dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
BM.dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(BM.dat),c("Lung_gene","Brain_gene","BMS_update","immunecheck","immune_act","immunsupp"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

# make a annotation 
ann_col <- mod[,c(6:8)]
ann_col <- data.frame(apply(ann_col,2,function(x){ifelse(x<quantile(x,0.33),"Low",ifelse(x>quantile(x,0.67),"High","Medium"))}))
# immunegene <- c("ADORA2A","ARHGEF5","BTLA","CD160","CD244","CD27","CD274","CD276","CD47","CD80","CEACAM1","CTLA4",
#     "GEM","HAVCR2","ICOS","IDO1","LAG3","PDCD1","TNFSF4","VISTA","VTCN1")
# gene.f <- immunegene[which(immunegene%in%rownames(BM.dat))]
gene.f <- c("LGALS9","ARG1","CD47","IDO1","VEGFA","IL10","TGFB1","PVR") # CD274 not found 
gene.f <- c("LGALS9","ARG1","CD47","IDO1","VEGFA","IL10","CD274","TGFB1","PVR") # CD274 should in E-MTAB 
dat <- BM.dat[gene.f,]


tmp.ann <- t(sapply(gene.f,function(x){
    c(cor.test(mod$Brain_gene_norm,as.numeric(dat[x,]),method="spearman")$p.value,
    cor.test(mod$Lung_gene_norm,as.numeric(dat[x,]),method="spearman")$p.value,
    cor.test(mod$BMS_update_norm,as.numeric(dat[x,]),method="spearman")$p.value
    )
}))
ann_row <- data.frame(apply(tmp.ann,2,function(x){ifelse(x>0,"pos","neg")}))[,,drop=F]
colnames(ann_row) <- c("Brain.cor","Lung.cor","Stem.cor")

ann_colors = list(
    Brain.cor = c(pos = "steelblue", neg = "red"),
    Lung.cor = c(pos = "steelblue", neg = "red"),
    Stem.cor = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = colorRampPalette(c('white',"#44AC48"))(100),
    Brain_gene_norm = colorRampPalette(c('white',"#269DCB"))(100),
    BMS_update_norm = colorRampPalette(c('white',"#F6BB3D"))(100)
)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_immunesuppression_gene_heatmap.pdf",useDingbats=F,width=15,height=12)
pheatmap(dat[],annotation_col=ann_col,scale="row",breaks=seq(-5,5,0.1),annotation_row=ann_row,
    annotation_colors=ann_colors,color=colorRampPalette(c('steelblue',"steelblue",'white',"red",'red'))(100),
    cellwidth=8,cellheight=12)
dev.off()











#=====================================================================================================================================================
# 2021-8-17
# MSGDB GO BP immune response activation signature 
# do this signature with Lung_like 
#=====================================================================================================================================================
# this signature seems not good 
# tmp <- as.vector(read.table("/public/workspace/lily/tmp/GOBP_immune_act.txt")[,1])
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(tmp,"immune_act",out="/public/workspace/lily/MOD_file/immune_act.mod") # make a mod file 
# calculate signature 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","Brain_gene","BMS_update","immunecheck","immune_act"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)


# E-MTAB
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","Brain_gene","BMS_update","immunecheck","immune_act"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)







#===================================================================================================================================================
# 2021-8-18
# 2021-10-14 2021-10-15
# use GSE14108 to plot heatmap
# T cell marker and Myeloid marker
#===================================================================================================================================================
library(pheatmap)
Tgene <- c("CD3D","CD3E","CD2","CD96")
Mgene <- c("ITGAM","CD14","CD163","CCR2")


ann_row <- mod[,6:8,drop=F]
ann_row <- data.frame(apply(ann_row,2,function(x){ifelse(x<quantile(x,0.33),"Low",ifelse(x>quantile(x,0.67),"High","Medium"))}))
#ann_row <- ann_row[order(ann_row$Lung_gene_norm),,drop=F]

tmp.dat <- dat[c(Tgene,Mgene),]
tmp.dat <- t(tmp.dat)

dat.BM <- dat
gene <- c(Tgene,Mgene)
dat <- dat.BM[gene,]
tmp.ann <- t(sapply(gene,function(x){
    c(cor.test(mod$Brain_gene_norm,as.numeric(dat[x,]),method="spearman")$estimate,
    cor.test(mod$Lung_gene_norm,as.numeric(dat[x,]),method="spearman")$estimate,
    cor.test(mod$BMS_update_norm,as.numeric(dat[x,]),method="spearman")$estimate
    )
}))
ann_col <- data.frame(apply(tmp.ann,2,function(x){ifelse(x>0,"pos","neg")}))[,,drop=F]
colnames(ann_col) <- c("Brain.cor","Lung.cor","Stem.cor")
ann_colors = list(
    Brain.cor = c(pos = "steelblue", neg = "red"),
    Lung.cor = c(pos = "steelblue", neg = "red"),
    Stem.cor = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = colorRampPalette(c('white',"#44AC48"))(100),
    Brain_gene_norm = colorRampPalette(c('white',"#269DCB"))(100),
    BMS_update_norm = colorRampPalette(c('white',"#F6BB3D"))(100)
)


ann_colors = list(
    Brain = c(pos = "steelblue", neg = "red"),
    Lung = c(pos = "steelblue", neg = "red"),
    Stem = c(pos = "steelblue", neg = "red"),
    Lung_gene_norm = c(High = "red", Medium="white",Low = "green"),
    Brain_gene_norm = c(High = "red", Medium="white",Low = "green"),
    BMS_update_norm = c(High = "red", Medium="white",Low = "green")
)






pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_myeloid_Tcell_heatmap.pdf",useDingbats=F)
pheatmap(tmp.dat,annotation_row=ann_row,scale="column",
    breaks=seq(-3,2,0.1),cellwidth=8,cellheight=8,annotation_col=ann_col,annotation_colors=ann_colors,
    color=colorRampPalette(c('steelblue',"steelblue",'white',"red",'red'))(60))
dev.off()


pheatmap(tmp.dat[,rownames(ann_col)],annotation_col=ann_col,scale="row",cluster_cols=F,cluster_rows=F,breaks=seq(-3,2,0.1),
    color=colorRampPalette(c('green',"green",'black',"red",'red'))(60))
















#==============================================================================================================================================
# 2021-8-20
# plot immune activation cytokine result 
# not siginificant 
#=============================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Lung_gene","Brain_gene","BMS_update","Treg","Cytokine_act","immunecheck"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

# plot result in Supplementary 
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE14108_Lung.Cytokine.pdf",useDingbats=F)
plot(mod$Lung_gene_norm,mod$Cytokine_act_norm,main="Lung like tumor with act Cytokine",xlab="Lung cancer signature",ylab="act Cytokine signature",pch=19)
dev.off()














#=============================================================================================================================================
# 2021-10-5
# LCBM check immune supressive gene expression percentage 
#=============================================================================================================================================
library(Seurat)
gene <- c("CD274","CD47","IDO1","LGALS9","IL10","TGFB1","VEGFA","ARG1","PVR")

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$group)

percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
           res.list[[i]] <- c(0,0,0,0)
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}

res.dat.list <- percent_feature(dat,gene,group="group")

list2mat <- function(res.list){
    res.dat <- matrix(unlist(res.list),ncol=length(res.list[[1]]),byrow=T)
    rownames(res.dat) <- unique(sapply(strsplit(as.vector(names(unlist(res.list))),"\\."),function(x){x[[1]]}))
    colnames(res.dat) <- unique(sapply(strsplit(as.vector(names(unlist(res.list))),"\\."),function(x){x[[2]]}))
    return(res.dat)
}

mat <- list2mat(res.dat.list)

# plot result 
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_subtumor_immunesupressive_percent.pdf",useDingbats=F)
pheatmap::pheatmap(mat,scale="row",show_colnames=T,color=colorRampPalette(c("steelblue",'white',"#E41A1C"))(100),cellwidth=10,cellheight=10)
dev.off()

# by the way plot other samples 
library(pheatmap)
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/LCBM_all_immunesuppressive_percent_heatmap.pdf",useDingbats=F)
pheatmap(mat.e[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.e")
pheatmap(mat.d[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.d")
# pheatmap(mat.lx701[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
# pheatmap(mat.lx681[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
# pheatmap(mat.lx255b[,1:3],scale="row",cluster_rows=F,cluster_cols=F)
pheatmap(mat.lcbm[,c(3,1,2)],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.lcbm")
pheatmap(mat.gse123902[,1:3],scale="row",cluster_rows=F,cluster_cols=F,cellwidth=10,cellheight=10,color=colorRampPalette(c('steelblue','white',"#E41A1C"))(100),main="mat.gse123902")
dev.off()

















#============================================================================================================================================
# 2021-11-23
# percentage of BMS.group tumor for 4 samples
#============================================================================================================================================

percent_feature <- function(dat,genelist,group){
    res.list <- c()
    for(i in 1:length(genelist)){
        dat$tmp_gene <- ifelse(dat[["RNA"]]@data[genelist[i],]>0,"Y","N")
        if(all(dat$tmp_gene=="N")){
           res.list[[i]] <- rep(0,length=length(table(dat@meta.data[,group])))
        }else{
           res.list[[i]] <- apply(table(dat$tmp_gene,dat@meta.data[,group]),2,function(x){x/sum(x)})[2,] 
        }
        
        names(res.list)[i] <- genelist[i]
    }
    return(res.list)
}

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
tumor.d.per <- percent_feature(dat,c("CCL2","CSF1"),"BMS.group")
dat@active.ident <- factor(dat$BMS.group)
tumor.d.exp <- AverageExpression(dat,assays="RNA",features=c("CCL2","CSF1"))$RNA

tumor.e <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
tumor.e.per <- percent_feature(tumor.e,c("CCL2","CSF1"),"BMS.group")
tumor.e@active.ident <- factor(tumor.e$BMS.group)
tumor.e.exp <- AverageExpression(tumor.e,assays="RNA",features=c("CCL2","CSF1"))$RNA

gse123902 <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
DefaultAssay(gse123902) <- "RNA"
gse123902.per <- percent_feature(gse123902,c("CCL2","CSF1"),"BMS.group")
gse123902@active.ident <- factor(gse123902$BMS.group)
gse123902.exp <- AverageExpression(gse123902,assays="RNA",features=c("CCL2","CSF1"))$RNA

lcbm <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
DefaultAssay(lcbm) <- "RNA"
lcbm.per <- percent_feature(lcbm,c("CCL2","CSF1"),"BMS.group")
lcbm@active.ident <- factor(lcbm$BMS.group)
lcbm.exp <- AverageExpression(lcbm,assays="RNA",features=c("CCL2","CSF1"))$RNA


tmp.res <- data.frame(
sample = rep(c("tumor.e","tumor.d","gse123902","lcbm"),each=2),
CCL2.exp = c(as.numeric(tumor.e.exp["CCL2",c(1,3)]),as.numeric(tumor.d.exp["CCL2",c(1,3)]),as.numeric(gse123902.exp["CCL2",c(1,3)]),as.numeric(lcbm.exp["CCL2",c(1,3)])),
CCL2.per = c(as.numeric(tumor.e.per$CCL2[c(1,3)]),as.numeric(tumor.d.per$CCL2[c(1,3)]),as.numeric(gse123902.per$CCL2[c(1,3)]),as.numeric(lcbm.per$CCL2[c(1,3)])),
CSF1.exp = c(as.numeric(tumor.e.exp["CSF1",c(1,3)]),as.numeric(tumor.d.exp["CSF1",c(1,3)]),as.numeric(gse123902.exp["CSF1",c(1,3)]),as.numeric(lcbm.exp["CSF1",c(1,3)])),
CSF1.per = c(as.numeric(tumor.e.per$CSF1[c(1,3)]),as.numeric(tumor.d.per$CSF1[c(1,3)]),as.numeric(gse123902.per$CSF1[c(1,3)]),as.numeric(lcbm.per$CSF1[c(1,3)])),
group = rep(c("high","low"),times=4)
)

tmp.res <- reshape2::melt(tmp.res)
res <- data.table::rbindlist(
    list( cbind(tmp.res[1:8,],tmp.res[9:16,4]) , cbind(tmp.res[17:24,],tmp.res[25:32,4]) ),
    use.names=FALSE 
)
colnames(res) <- c("sample","group","gene","exp","percent")
res$gene <- sapply(strsplit(as.vector(res$gene),"\\."),function(x){x[[1]]})
res <- as.data.frame(res,stringsAsFactors=F)



# plot result
library(ggplot2)
res$class <- paste0(res$sample,"-",res$gene)

tmp.percent <- c()
for(i in c(1,3,5,7,9,11,13,15)){
   tmp.percent <- c(tmp.percent,as.numeric(scale(res$percent[c(i,i+1)],center=F)))
}
res$scale.percent <- tmp.percent

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/CCL2_CSF1_BMS_high_low.pdf",useDingbats=F)
ggplot(res,aes(x=group,y=(class))) + geom_point(aes(size=scale.percent,color=exp)) + scale_colour_gradientn(colours=c("steelblue","#F1AB3E"))
dev.off()






#============================================================================================================================================
# 2021-11-23
# GSE14108 use CCL2 and CSF1 and MDM-MG
#============================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")
# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","BMDM_marker","MG_marker"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$MDM.MG <- mod$BMDM_marker_norm - mod$MG_marker_norm 
mod$CCL2 <- as.numeric(dat["CCL2",])
mod$CSF1 <- as.numeric(dat["CSF1",])

cor.test(mod$CCL2,mod$BMS_update_norm,method="spearman")
cor.test(mod$CCL2,mod$BMDM_marker_norm,method="spearman")
cor.test(mod$CCL2,mod$MG_marker_norm,method="spearman")
cor.test(mod$CCL2,mod$MDM.MG,method="spearman")
cor.test(mod$CSF1,mod$MDM.MG,method="spearman")




# plot result
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/CCL2_BMS_GSE14108.pdf",useDingbats=F)
par(mfrow=c(2,2))
plot(mod$BMS_update_norm,mod$CCL2,main="CCL2 with BMS",pch = 20, cex = 2.5)
abline(lm(mod$CCL2~mod$BMS_update_norm),col="red")

plot(mod$BMDM_marker_norm,mod$CCL2,main="CCL2 with BMDM",pch = 20, cex = 2.5)
abline(lm(mod$CCL2~mod$BMDM_marker_norm),col="red")

plot(mod$MG_marker_norm,mod$CCL2,main="CCL2 with MG",pch = 20, cex = 2.5)
abline(lm(mod$CCL2~mod$MG_marker_norm),col="red")

plot(mod$MDM.MG,mod$CCL2,main="CCL2 with MDM.MG",pch = 20, cex = 2.5)
abline(lm(mod$CCL2~mod$MDM.MG),col="red")
dev.off()

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig5/CSF1_BMS_GSE14108.pdf",useDingbats=F)
par(mfrow=c(2,2))
plot(mod$BMS_update_norm,mod$CSF1,main="CSF1 with BMS",pch = 20, cex = 2.5)
abline(lm(mod$CSF1~mod$BMS_update_norm),col="red")

plot(mod$BMDM_marker_norm,mod$CSF1,main="CSF1 with BMDM",pch = 20, cex = 2.5)
abline(lm(mod$CSF1~mod$BMDM_marker_norm),col="red")

plot(mod$MG_marker_norm,mod$CSF1,main="CSF1 with MG",pch = 20, cex = 2.5)
abline(lm(mod$CSF1~mod$MG_marker_norm),col="red")

plot(mod$MDM.MG,mod$CSF1,main="CSF1 with MDM.MG",pch = 20, cex = 2.5)
abline(lm(mod$CSF1~mod$MDM.MG),col="red")
dev.off()
































































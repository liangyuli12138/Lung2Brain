

# 2022-5-23 
# prepare for some public data 
# GSE200563
# paired LCBM and MLUAD data 
library(Biobase)
tmp <- GEOquery::getGEO(filename="~/metastasis/data/verify/GSE200563/GSE200563_series_matrix.txt.gz",getGPL=F) 
sampleinfo <- pData(tmp)[,c(8,26,41:47)]
colnames(sampleinfo) <- c("group","description","metastasis_diagnosis_to_death_m",
    "gender","grade","histological","location_of_brm","metastasis_intervals_to_brain",
    "primary_diagnosis_death_m"
    )

# recode 
sampleinfo$group <-  car::recode(sampleinfo$group,
	" 'Metastatic lung cancer in a lymph node'='MLN';
	'Metastatic lung cancer in the brain'='BM';
	'non-tumor brain'='NormalB';
	'primary lung cancer'='MLUNG';
    'tumor brain microenvironment in metastatic brain'='TME_BM';
	'tumor immune microenvironment in the metastatic brain'='TIME_BM';
	'tumor immune microenvironment in the primary lung cancer'='TIME_LUNG' "
	)

snum <- read.table("~/metastasis/data/verify/GSE200563/GSE200563_sample",sep="\t")
all(snum$V1==rownames(sampleinfo))
sampleinfo$Patient <- sapply(strsplit(as.vector(snum[,2]),"\\s"),function(x){paste(x[1:2],collapse="_")})
sampleinfo$Patient[1:7] <- paste("Nt_B_",1:7)
saveRDS(sampleinfo,file="~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# get data 
data <- read.table("~/metastasis/data/verify/GSE200563/GSE200563_series_matrix.txt.gz",comment.char="!",sep="\t",header=T)
rownames(data) <- data$ID_REF
data$ID_REF <- NULL
saveRDS(data,file="~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")






# make a combination for TCGA and this sample
# TCGA just use early stage sample 
a <- readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD.clin.f.RDS")
tmp.tcga.luad <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
tcga.luad <- tmp.tcga.luad[,a$Row.names[which(a$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage II","Stage IIA","Stage IIB"))]]
# saveRDS(tcga.luad,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp_StageI_II.RDS")
# filter some sample 
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# need filter 
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM","MLUNG")),]
dat <- tmp.dat[,rownames(sampleinfo)]

merge.dat <- merge(tcga.luad,dat,by="row.names")
rownames(merge.dat) <- merge.dat[,1]
merge.dat$Row.names <- NULL

res.tcga <- merge.dat[,1:442]
res.gse <- merge.dat[,443:ncol(merge.dat)]
res.gseinfo <- sampleinfo

res <- list(TCGA_Data=res.tcga,GSE_Data=res.gse,GSE_info=res.gseinfo)
saveRDS(res,file="~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")
# 2022-5-24 have try combat and house-keeping gene not ok
#######################################################################################################################
# grep("PGK1",rownames(merge.dat))
# res.dat <- apply(merge.dat,2,function(x){x/x[10195]})

# batch <- data.frame(row.names=colnames(merge.dat),sample=colnames(merge.dat),batchinfo=rep("TCGA",ncol(merge.dat)),stringsAsFactors=F)
# batch$batchinfo[which(rownames(batch)%in% rownames(sampleinfo))] <- "GSE"

# edata <- sva::ComBat(dat=as.matrix(merge.dat), batch=batch$batchinfo)


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_V6","chr5p12","chr5p13","chr5p14","chr5p15"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$group <- "TCGA"
mod$group[which(rownames(mod)%in%rownames(sampleinfo)[which(sampleinfo$group=="BM")])] <- "BM"
mod$group[which(rownames(mod)%in%rownames(sampleinfo)[which(sampleinfo$group=="MLUNG")])] <- "MLUNG"



#============================================================================================================================
# GSE198291
library(Biobase)
tmp <- GEOquery::getGEO(filename="~/metastasis/data/verify/GSE198291/GSE198291_series_matrix.txt.gz",getGPL=F) 
sampleinfo <- pData(tmp)[,c(8,39)]
colnames(sampleinfo) <- c("source","type")
sampleinfo$Patient <- sapply(strsplit(as.vector(sampleinfo$source),"_"),function(x){x[[1]]})













#============================================================================================================================
#============================================================================================================================
# also change GSE14108 data with TCGA LUAD early stage data 
tcga.luad <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp_StageI_II.RDS")
# GSE14108 data 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")

merge.dat <- merge(tcga.luad,dat,by="row.names")
rownames(merge.dat) <- merge.dat[,1]
merge.dat$Row.names <- NULL

res.tcga <- merge.dat[,1:442]
res.gse <- merge.dat[,443:ncol(merge.dat)]
res <- list(TCGA_Data=res.tcga,GSE_Data=res.gse)

saveRDS(res,file="~/metastasis/data/verify/GSE14108/TCGA_stage_I_II_GSE14108_BM_MLUAD_list.RDS")










#============================================================================================================================
#============================================================================================================================
# Prepare for MSK cell 2021 metastasis panel
# 2022-5-25
tmp.seg <- read.table("/public/workspace/lily/metastasis/data/MSK_Cell_2021/msk_met_2021_segments.seg",header=T,sep="\t")
dat <- tmp.seg[which(tmp.seg$chrom==5&tmp.seg$loc.start<46100000&tmp.seg$loc.end<46100000),]
dat$length <- dat$loc.end - dat$loc.start
dat$value <- dat$seg.mean * dat$length
tmp.res <- aggregate(.~ID,data=dat[,c(1,7,8)],FUN=sum)
tmp.res$res <- tmp.res$value / tmp.res$length

# load sampleinfo 
tmp.info <- read.table("/public/workspace/lily/metastasis/data/MSK_Cell_2021/msk_met_2021_clinical_data.tsv",sep="\t",header=T)
idx.1 <- tmp.info$Sample.ID[which(tmp.info$Metastatic.patient=="FALSE"&tmp.info$Sample.Type=="Primary")] #原发未转
idx.2 <- tmp.info$Sample.ID[which(tmp.info$Metastatic.Site=="CNS/Brain"&tmp.info$Sample.Type=="Metastasis")] #转移灶
idx.3 <- tmp.info$Sample.ID[which(tmp.info$Distant.Mets..CNS.Brain=="Yes"&tmp.info$Sample.Type=="Primary")] # 转移原发

median(tmp.res$res[which(tmp.res$ID %in% idx.1)])
median(tmp.res$res[which(tmp.res$ID %in% idx.2)])
median(tmp.res$res[which(tmp.res$ID %in% idx.2)])


# plot a result 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/MSK_Cell_2021_LUAD_Brain_seg.pdf",useDingbats=F)
boxplot(tmp.res$res[which(tmp.res$ID %in% idx.1)],tmp.res$res[which(tmp.res$ID %in% idx.3)],
tmp.res$res[which(tmp.res$ID %in% idx.2)],
name=c("nMLUAD","MLUAD","BM"),
outline=F) #170,111,207
dev.off()



# wilcox.test(tmp.res$res[which(tmp.res$ID %in% idx.1)],tmp.res$res[which(tmp.res$ID %in% idx.2)])



















# 2022-8-16
# GSE166720
# not OK
#================================================================================================================================================
library(Biobase)
tmp <- GEOquery::getGEO(GEO="GSE166720",getGPL=F) 
sampleinfo <- pData(tmp[[1]])[,c(1,49:51)]
colnames(sampleinfo) <- c("title","subtype","Tstage","tissue")
sampleinfo$subtype[which(sampleinfo$subtype=="indolent")] <- "Indolent"
saveRDS(sampleinfo,file="/public/workspace/lily/metastasis/data/verify/GSE166720/GSE166720_sampleinfo.RDS")

tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE166720/GSE166720_NIHAD_FPKM_normalized.txt",sep="\t",header=T)
rownames(tmp.dat) <- tmp.dat$X
tmp.dat$X <- NULL
all(colnames(tmp.dat)==sampleinfo$title)
sampleinfo$OSMR <- as.numeric(tmp.dat["OSMR",])
sampleinfo$OSM <- as.numeric(tmp.dat["OSM",])





# 2022-8-16
# GSE162698
# not significant but OK
#================================================================================================================================================

tmp.dat <- data.table::fread("~/metastasis/data/verify/GSE162698/GSE162698_tpm_summary_invitro_Macs.txt")
tmp.dat[6,3:17] <- tmp.dat[5,3:17]
dat <- tmp.dat[6:nrow(tmp.dat),]
colnames(dat) <- as.character(dat[1,])
dat <- dat[-1,]
dat <- as.data.frame(dat)
dat$ID <- NULL

dat.t <- data.frame(apply(dat[,-1],2,function(x){as.numeric(as.character(x))}))
dat.t$symbol <- dat$symbol

tmp.res <- aggregate(.~symbol,data=dat.t,FUN=max)
rownames(tmp.res)<- tmp.res$symbol
tmp.res <- tmp.res[-1,]

saveRDS(tmp.res,file="~/metastasis/data/verify/GSE162698/GSE162698_tpm.RDS")







# GSE162669
# not significant but OK
#================================================================================================================================================

tmp.dat <- data.table::fread("~/metastasis/data/verify/GSE162669/GSE162669_tpm_summary_exvivo_AMTAM.txt")
tmp.dat[5,3:17] <- tmp.dat[3,3:17]
dat <- tmp.dat[5:nrow(tmp.dat),]
colnames(dat) <- as.character(dat[1,])
dat <- dat[-1,]
dat <- as.data.frame(dat)
dat$ID <- NULL

dat.t <- data.frame(apply(dat[,-1],2,function(x){as.numeric(as.character(x))}))
dat.t$symbol <- dat$symbol

tmp.res <- aggregate(.~symbol,data=dat.t,FUN=max)
rownames(tmp.res)<- tmp.res$symbol
tmp.res <- tmp.res[-1,]
tmp.res$symbol <- NULL

saveRDS(tmp.res,file="~/metastasis/data/verify/GSE162669/GSE162669_tpm.RDS")



# 2022-8-16
# GSE116946
# not significant 
#================================================================================================================================================

tmp.dat <- read.csv("~/metastasis/data/verify/GSE116946/GSE116946_NormCounts_DESeq_EGM_UoS.csv")
rownames(tmp.dat) <- tmp.dat$Gene

	samplenum = 74
    gene.info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",header=T,sep="\t")
    gene.info$length <- gene.info$Gene.end..bp. -gene.info$Gene.start..bp.
    colnames(gene.info)[1] <- "gene_ensemble"

    res.tmp <- merge(tmp.dat,gene.info,by.x="Gene",by.y="Gene.name")
    rownames(res.tmp) <- res.tmp$gene_ensemble
    res.tmp[,c(2:(samplenum+1),(samplenum+5))] -> mat
    mat.tmp <- apply(mat,1,function(x){x/x[length(x)]})
    mat.res <- apply(mat.tmp,1,function(x){(x/sum(x))*10^6})
    mat.res <- data.frame(mat.res)
    mat.res$gene_ensemble <- rownames(mat.res)
    res <- merge(mat.res,gene.info[,c(1,2)],by="gene_ensemble")
    res$gene_ensemble <- NULL
    res$length <- NULL

    aggregate(.~Gene.name,data=res,FUN=median) -> res.final
    rownames(res.final) <- res.final$Gene.name
    res.final$Gene.name <- NULL
    #head(res.final)




dat <- res.final
sampleinfo <- data.frame(row.names=colnames(dat),Sample=colnames(dat))
sampleinfo$Group <- gsub("\\.[0-9]+$","",sampleinfo$Sample)
sampleinfo$SID <- sapply(strsplit(as.character(sampleinfo$Sample),"\\."),function(x){x[length(x)]})



saveRDS(dat,file="~/metastasis/data/verify/GSE116946/GSE116946_dat.RDS")
saveRDS(sampleinfo,file="~/metastasis/data/verify/GSE116946/GSE116946_sampleinfo.RDS")









# GSE100412
# Tumor cells 
#================================================================================================================================================
tmp.dat <- read.csv("~/metastasis/data/verify/GSE100412/GSE100412_CMTandLLC_FPKM_table_V2.txt",sep="\t",header=T)
colnames(tmp.dat) <- c("Gene","location","LLC_InVitro_1","LLC_InVitro_2","LLC_InVitro_3","LLC_InVitro_4","LLC_InVitro_5",
"LLC_InVivo_1","LLC_InVivo_2","LLC_InVivo_3","LLC_InVivo_4","LLC_InVivo_5",
"LLC_Pass_1","LLC_Pass_2","LLC_Pass_3","LLC_Pass_4","LLC_Pass_5",
"CMT_InVitro_1","CMT_InVitro_2","CMT_InVitro_3","CMT_InVivo_1","CMT_InVivo_2","CMT_InVivo_3","CMT_Pass_1","CMT_Pass_2","CMT_Pass_3")

tmp.dat$location <- NULL
rownames(tmp.dat) <- tmp.dat$Gene
tmp.dat$Gene <- NULL

sampleinfo <- data.frame(row.names=colnames(tmp.dat),Name=colnames(tmp.dat))
sampleinfo$Group <- sapply(strsplit(as.character(sampleinfo$Name),"_"),function(x){x[2]})
sampleinfo$Sample <- sapply(strsplit(as.character(sampleinfo$Name),"_"),function(x){paste0(x[1],"_",x[3])})


saveRDS(tmp.dat,file="~/metastasis/data/verify/GSE100412/GSE100412_FPKM.RDS")
saveRDS(sampleinfo,file="~/metastasis/data/verify/GSE100412/GSE100412_sampleinfo.RDS")















# 2022-9-23
# analysis about organsim
# GSE8549
#==================================================================================================================================================
load("~/metastasis/data/verify/GSE18549/ALL/data.RData")
rownames(type)<- type$GSM_num
all(rownames(type)==colnames(dat))

# 1. Lung to 
# grep("Lung",type$cancer,value=T)
# lung to kidney just 1 sample
dat.f <- dat[,grep("Lung",type$cancer)[c(8,10:17)]]
type.f <- type[grep("Lung",type$cancer)[c(8,10:17)],]

# calculate BMS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$type <- type.f$cancer
mod$OSMR <- as.numeric(dat.f["OSMR",])


# 2. all metastasis
all(rownames(type)==colnames(dat))
dat.f <- dat[,which(type$cancer%in%names(which(table(type$cancer)>2)))]
type.f <- type[which(type$cancer%in%names(which(table(type$cancer)>2))),]

# calculate BMS 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.f),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)

mod$type <- type.f$cancer
mod$OSMR <- as.numeric(dat.f["OSMR",])

library(ggplot2)
library(forcats)
library(dplyr)
tmp.res <- aggregate(.~type,data=mod,FUN=median)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE18549_OSMR.pdf",useDingbats=F)
mod %>% mutate(type = fct_reorder(type,OSMR,.fun="median")) %>%
ggplot(aes(x=type,y=OSMR,fill=type))+geom_boxplot() + theme_classic()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
dev.off()

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE18549_BMS.pdf",useDingbats=F)
mod %>% mutate(type = fct_reorder(type,BMS_update_norm,.fun="median")) %>%
ggplot(aes(x=type,y=BMS_update_norm,fill=type))+geom_boxplot() + theme_classic()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))
dev.off()





# 2022-9-23
#=====================================================================================================================================================
# GSE123902
filepath <- grep("^GSM[0-9]+_MSK",dir("~/metastasis/data/verify/GSE123904"),value=T)
tmp.list <- list()
for(i in 1:length(filepath)){
    tmp.dat <- read.csv(file=paste0("~/metastasis/data/verify/GSE123904/",filepath[i]))
    samplename <- paste(strsplit(as.character(filepath[i]),"_")[[1]][3:4],collapse="_")
    rownames(tmp.dat) <- paste0(samplename,"_",tmp.dat$X)
    tmp.dat$X <- NULL
    

    library(Seurat)
    tmp_dat <- CreateSeuratObject(count=t(tmp.dat),project=samplename,assay="RNA")
    tmp_dat = NormalizeData(object = tmp_dat)
    tmp_dat <- FindVariableFeatures(object = tmp_dat)
    all.genes <- rownames(x = tmp_dat)
    tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)

    tmp.list[[i]] <- tmp_dat

}

integration.anchors <- FindIntegrationAnchors(object.list = tmp.list)
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

saveRDS(inte,file="~/metastasis/data/verify/GSE123904/GSE123904_17S_human_dat.RDS")
DefaultAssay(inte) <- "RNA"
FeaturePlot(inte,features=c("EPCAM","KRT19","CDH1"),label=T,label.size=6)

inte$group <- sapply(strsplit(colnames(inte),"_"),function(x){paste0(x[1],"_",x[2])})
inte$type <- sapply(strsplit(colnames(inte),"_"),function(x){x[2]})

inte$type[which(inte$orig.ident%in%c("LX255B","LX681","LX701"))] <- "Lung2Brian"
inte$type[which(inte$orig.ident%in%c("LX699"))] <- "Lung2Adrenal"
inte$type[which(inte$orig.ident%in%c("LX666"))] <- "Lung2Bone"
saveRDS(inte,file="~/metastasis/data/verify/GSE123904/GSE123904_17S_human_dat.RDS")
# 11 ,12
#=======================================================================================================================================
library(Seurat)
tmp.dat <- readRDS("~/metastasis/data/verify/GSE123904/GSE123904_17S_human_dat.RDS")
DefaultAssay(tmp.dat) <- "RNA"

dat <- subset(tmp.dat,cells=which(tmp.dat$seurat_clusters%in%c(11,12)))




dat$OSMR <- ifelse(dat[["RNA"]]@data["OSMR",]>0,"Y","N")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)










#=======================================================================================================================
# 2022-10-5
# GSE198291
# analysis scRNA 
library(Seurat)
tmp.dat <- read.csv("~/metastasis/data/verify/GSE198291/GSE198291_allcounts2.csv")
rownames(tmp.dat) <- tmp.dat$X
tmp.dat$X <- NULL

# make seurat object
dat <- CreateSeuratObject(counts=tmp.dat,assay="RNA")
# seurat object
	tmp_dat = NormalizeData(object = dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.8)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
dat <- tmp_dat
dat$group <- substr(sapply(strsplit(colnames(dat),"_"),function(x){x[2]}),1,3)
dat$sample <- sapply(strsplit(colnames(dat),"_"),function(x){x[1]})
# saveRDS(dat,file="~/metastasis/data/verify/GSE198291/GSE198291_dat.RDS")

FeaturePlot(dat,features=c("PTPRC","EPCAM","CDH1","CD68"),label=T,label.size=6)
# EPCAM 4,6,7,8

dat <- readRDS("~/metastasis/data/verify/GSE198291/GSE198291_dat.RDS")
subdat <- subset(dat,cells=which(dat$seurat_clusters%in%c(4,6,7,8)))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("BMS_update"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)











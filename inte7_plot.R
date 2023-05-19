
#================================================================================================================
# plot for result and report 
#================================================================================================================
# signature is BMS_test ,including 57 genes 
#================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS") 

#          B_cell     Endothelial      Fibroblast      maliganant         Myeloid
#            1884             392            1046           10321            9429
# Oligodendrocyte          T_cell          unknow
#             954           11874             180
cols <- c('#377EB8','#910241','#984EA3','#E41A1C','#F29403','#8CA77B','#B2DF8A','#999999')

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/lanscape.pdf")
DefaultAssay(dat) <- "RNA"
DimPlot(dat,label=F)
DimPlot(dat,group.by="type",cols=cols)
FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=F,order=F) # T cell
FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','FCGR1A'),label=F,order=F) # Myeloid
FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=F,order=T) # B cell 
FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=F,order=T) # Oligodendrocyte
FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=F,order=T) # Fibroblast/Vascular
FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=F,order=T) # Endothelial
FeaturePlot(dat,features=c("EGFR","EPCAM","PTPRC"),label=F,order=T,cols=c("lightgrey", "red"))
dev.off()

#==============================================================================================================
#
# genes 
#==============================================================================================================
data <- dat[['RNA']]@data

tumor_gene <- apply(data[,which(dat$maliganant=="tumor")],2,function(x){length(which(x[]>0))})
ntumor_gene <- apply(data[,which(dat$maliganant=="non-tumor")],2,function(x){length(which(x[]>0))})

pdf("./tmp.pdf")
boxplot(tumor_gene,ntumor_gene,names=c('malignant','non-malignant'),main="express gene numbers",col=c("#E41A1C","#A6CEE3"),outline=F)
legend('topright',legend=paste('p =',round(wilcox.test(tumor_gene,ntumor_gene)$p.value,5)),bty='n')
dev.off()











#============================================================================================================
# 
# BM tumor and LC tumor plot to sketch 
#============================================================================================================

BM <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_BM.RDS")
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/BM_sigtumor.pdf")
DimPlot(BM,cells.highlight=colnames(BM)[which(BM$seurat_clusters%in%c(6,9))])
dev.off()


LC <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LC_sigtumor.pdf")
DimPlot(LC,cells.highlight=colnames(LC)[which(LC$seurat_clusters%in%c(0,11))])
dev.off()







#============================================================================================================
#
# cell cycle 
#============================================================================================================
dat <- CellCycleScoring(dat,s.features=cc.genes.updated.2019$s.genes,g2m.features=cc.genes.updated.2019$g2m.genes)
subdat <- subset(dat,cells=which(dat$maliganant=="tumor"))

apply(a,2,function(x){x/sum(x)})-> tmp

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/cellcyle.pdf")
barplot(tmp,col=c('#54B0E4','#1B9E77','#B2DF8A'))
legend("topleft",pch=c(15,15,15),legend=c("G1","G2M","S"),col=c('#54B0E4','#1B9E77','#B2DF8A'),bty="n",horiz=TRUE)
dev.off()







#================================================================================================================
# bar plot
# stack  bar 
#================================================================================================================
table(dat$type,dat$type_group) ->a
apply(a,2,function(x){x/sum(x)}) -> tmp
barplot(tmp,col=c('#377EB8','#910241','#984EA3','#E41A1C','#F29403','#8CA77B','#B2DF8A','#999999'))





#================================================================================================================
# TCGA stage 
#
#================================================================================================================
luad_clin <- read.table('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt',sep='\t',header=T)
luad_clin <- luad_clin[,c(1,102:105)]
rownames(luad_clin) <- gsub("-",".",luad_clin$sampleID)
load("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_data.RData") # dat 
dat -> luad_dat
#luad_dat.f <- luad_dat[,which(colnames(luad_dat)%in%rownames(luad_clin)[grep("Stage",luad_clin$pathologic_stage)])]
#============================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
luad_mod <- mod.analyze2(as.matrix(luad_dat),"BMS_test","/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
#rownames(luad_mod) <- gsub("\\.","-",rownames(luad_mod))
dat.f <- merge(luad_clin,luad_mod,by="row.names")


#====================================== purity 
pur <- read.table("/public/workspace/lily/metastasis/data/verify/LUAD_absolute.txt",sep="\t",header=T)
pur$sampleID <- substr(pur$Sample.ID,1,15)
dat.ff <- merge(dat.f,pur,by="sampleID")



dat.fff <- dat.ff[grep("Stage",dat.f$pathologic_stage),]
dat.fff$Stage <- "stage IV"
dat.fff$Stage[which(dat.fff$pathologic_stage%in%c("stage I","Stage IA","Stage IB"))] <- "stage I"
dat.fff$Stage[which(dat.fff$pathologic_stage%in%c("Stage II","Stage IIA","Stage IIB"))] <- "stage II"
dat.fff$Stage[which(dat.fff$pathologic_stage%in%c("Stage IIIA","Stage IIIB"))] <- "stage III"

res <- dat.fff[which(dat.fff$ABSOLUTE>0.3),]

aggregate(BMS_test_norm~Stage,data=res,FUN=median)




#================================================================================================================
# use bulk RNA-seq to analysis stage 
# and use it to analysis 
#================================================================================================================
library(Seurat)

#====================================== GSE19804 ================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE19804/GSE19804_series_matrix.txt",sep="\t",header=T,comment.char="!") 
library(idmap1,lib.loc="/public/workspace/wangzy/R/x86_64-pc-linux-gnu-library/3.6.0")
ids=getIDs('gpl570')
GPL570 <- ids[,c("probe_id","symbol")]
colnames(GPL570)[1] <- "ID_REF"
dat.f <- merge(GPL570,dat,by="ID_REF")
#dat.ff <- dat.f[-grep("///",dat.f$symbol),]
dat.f$ID_REF <- NULL

dat.ff <- aggregate(.~symbol,data=dat.f,FUN=mean)
rownames(dat.ff) <- dat.ff$symbol
dat.ff$symbol <- NULL
saveRDS(dat.ff,file="/public/workspace/lily/metastasis/data/verify/GSE19804/expr.RDS")

#================================================================================================================
# load stage info 
#===============================================================================================================
stage <- read.table("/public/workspace/lily/metastasis/data/verify/GSE19804/stage_info.txt",sep='\t',header=T)
rownames(stage) <- stage$ID


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE19804/expr.RDS")
mod <- mod.analyze2(as.matrix(dat),"BMS_test","/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- data.frame(mod)
res <- merge(mod,stage,by="row.names")

res$state <- "-"
res$state[grep("1",res$Stage)] <- "I"
res$state[grep("2",res$Stage)] <- "II"
res$state[grep("3",res$Stage)] <- "III"
res$state[grep("4",res$Stage)] <- "IV"


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE19804_stage.pdf")
boxplot(BMS_test_norm~state,data=res,outline=F)
dev.off()







gene <- c('BRI3','ASPH','FKBP5','CEBPD','JUND','SKAP2','SPTSSA','MET','HOPX','TOP2A','EREG',
	'C4orf48','BIRC3','UPP1','GLS','APLP2','CPD','DCBLD2','FAM3C','TK1','IRX2','ERRFI1','MBNL1',
	'DUSP6','MAL2','TIPARP','BIRC5','HIST1H4C','SFTA2','ITPKA','MGLL','MKI67','S100A16','NFKBIA',
	'UBE2C','SLC20A1','HPGD','MUC1','S100A4','CXCL2','DUSP1','CENPF','VIM','ADK','LGALS1','CEBPB',
	'PGLS','PTTG1','CCNB1','PHLDA2','G0S2','SFTPB','NFKBIZ','IL8','LCN2','CCL20','KRT6A')

# library(pheatmap)
# stage$state <- "-"
# stage$state[grep("1",stage$Stage)] <- "I"
# stage$state[grep("2",stage$Stage)] <- "II"
# stage$state[grep("3",stage$Stage)] <- "III"
# stage$state[grep("4",stage$Stage)] <- "IV"

# stage <- stage[order(stage$state),]
# pheatmap(dat[which(rownames(dat)%in%gene),rownames(stage)],cluster_cols=F,annotation_col=stage,scale="row",color= colorRampPalette(c("steelblue",'white','red'))(500))






#================================================================================================================
# RABIT
#================================================================================================================
load('~/res_cor_inte/MOD/luad_mod.RData') # luad_mod
rabit_luad <- read.table('~/metastasis/data/TCGA_RABIT/LUAD/RABIT_LUAD.HiSeq.V2',sep='\t',header=T) #80 TFs 
rownames(rabit_luad) <- gsub("-",".",rabit_luad[,1])
rabit_luad$X <- NULL		#prepare

# get co-sample 
co_luad <- merge(luad_mod,rabit_luad,by='row.names')
rownames(co_luad) <- co_luad$Row.names
co_luad$Row.names <- NULL
tf_luad <- data.frame(t(apply(co_luad[,-c(1:3)],2,function(x,y){tmp=cor.test(x,y);c(tmp$estimate,tmp$p.value,length(which(x!=0))/nrow(co_luad))},y=co_luad$BMS_test_norm)))
# col 1:3 is mod-result 
colnames(tf_luad) <- c('cor','pvalue','percent')
tf_luad$log2Pvalue <- -log2(tf_luad$pvalue) # -log2(pvalue)
# set color
tf_luad$type <- 'up_regulate'
tf_luad$type[which(tf_luad$cor<0)] <- 'down_regulate'
# set label
tf_luad$label <-''
tf_luad$label[which(tf_luad$log2Pvalue>10)] <- rownames(tf_luad[which(tf_luad$log2Pvalue>10),])

pdf("./tf_LUAD.pdf")
ggplot(data=tf_luad,aes(x=cor,y=log2Pvalue,label=label,col=type))+geom_point()+geom_text(aes(label=tf_luad$label),check_overlap = TRUE,vjust = 0, nudge_y = 0.5, size=3.5)+theme_classic()+
	scale_colour_manual(values=c('steelblue1','red'))
dev.off()


#================================================================================================================
# GDSC 
#
dat <- read.table("/public/workspace/lily/Lung2Brain/inte7/GDSC/Cell_line_RMA_proc_basalExp.txt",sep="\t",header=T)
luad_sample <- read.table("/public/workspace/lily/Lung2Brain/inte7/GDSC/cellline_LUAD.txt",sep="\t",header=T)
tmp_dat <- dat[,c(1,which(colnames(dat)%in%paste0("DATA.",luad_sample$COSMIC_ID)))]
dat.f <- aggregate(.~GENE_SYMBOLS,data=tmp_dat,FUN=mean)
rownames(dat.f) <- dat.f$GENE_SYMBOLS
dat.f$GENE_SYMBOLS <- NULL
dat.f <- dat.f[-1,]
saveRDS(dat.f,file="/public/workspace/lily/Lung2Brain/inte7/GDSC/cellline169_exp.RDS")
#============================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
GDSC_mod <- mod.analyze2(as.matrix(dat.f),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)







#========================================================================================================================================
load("/public/workspace/lily/Lung2Brain/inte7/GDSC_mod.RData")
GDSC_mod <- data.frame(GDSC_mod)
GDSC_mod$COSMIC_ID <- gsub("DATA\\.","",rownames(GDSC_mod))
dat2 <- read.table("/public/workspace/lily/Lung2Brain/inte7/GDSC/GDSC2.txt",sep="\t",header=T)
dat2.f <- merge(dat2,luad_sample,by="COSMIC_ID")
dat2.ff <- merge(dat2.f,GDSC_mod,by="COSMIC_ID")
#========================================================================================================================================
res <- matrix(ncol=3)
drug <- as.vector(unique(dat2.ff$DRUG_NAME))
for(i in 1:length(drug)){
	tmp <- dat2.ff[which(dat2.ff$DRUG_NAME==drug[i]),]
	tmp1 <-cor.test(tmp$LN_IC50,tmp$BMS_test_norm)
	res <- rbind(res,c(tmp1$estimate,tmp1$p.value,drug[i]))
	#rownames(res)[nrow(res)] <- drug[i]
}

res <- data.frame(res)
res <- res[-1,]
colnames(res) <- c("Cor","p_value","drug_name")
res$Cor <- as.numeric(as.vector(res$Cor))
res$p_value <- as.numeric(as.vector(res$p_value))
res$p.adj <- p.adjust(res$p_value)

tmp <- res[which(res$p.adj<0.01),]
tmp$log2Pvalue <- (-log2(tmp$p.adj))
tmp$type <- "resistance"
tmp$type[which(tmp$Cor<0)] <- "sensitive"
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LUAD_drug.pdf")
ggplot(data=tmp,aes(x=Cor,y=log2Pvalue,label=drug_name,col=type))+geom_point()+geom_text(aes(label=tmp$drug_name),check_overlap = TRUE,vjust = 0, nudge_y = 0.5, size=3.5)+theme_classic()+
	scale_colour_manual(values=c('steelblue1','red'))
dev.off()









#================================================================================================================
# TCGA LUAD Mutation 
#
#================================================================================================================

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




#================================================================================================================
# mutation KRAS G12 point 
#
#================================================================================================================
dat <- read.table('~/metastasis/data/verify/TCGA_LUAD/mutation_broad_LUAD',sep='\t',header=T)
dat.f <- dat[which(dat$gene=="KRAS"),]
dat.ff <-dat.f[grep("p.G12",dat.f$Amino_Acid_Change),]
load("/public/workspace/lily/Lung2Brain/inte7/TCGA_LUAD_mod.RData")
luad_mod <- data.frame(luad_mod)
luad_mod$type <- "No_KRAS_mut"
luad_mod$type[which(rownames(luad_mod)%in%gsub("-",".",dat.f$sample))] <- "KRAS_mut"
luad_mod$type[which(rownames(luad_mod)%in%gsub("-",".",dat.ff$sample)&rownames(luad_mod)!="No_KRAS_G12_mut")] <- "KRAS_G12_mut"
aggregate(BMS_test_norm~type,data=luad_mod,FUN=median)





























#================================================================================================================
# run infercnv 
#===================================


library(infercnv)
dat <- readRDS("./Lung2Brain/inte7/inte7_ann.RDS")
subdat <- subset(dat,cells=which(dat$type%in%c("B_cell","maliganant")))

setwd('/public/workspace/lily/Lung2Brain/inte7/infercnv/')
dir.create('./B_cell_maliganat')
setwd('/public/workspace/lily/Lung2Brain/inte7/infercnv/B_cell_maliganat/')

write.table(subdat$type,paste('inte7',"cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
count<-as.matrix(subdat@assays$RNA@counts)
write.table(count,paste("inte7","count_exp.txt",sep='_'),sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste('inte7',"count_exp.txt",sep='_'),
         annotations_file=paste('inte7',"cell_info.txt",sep='_'),
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names=c("B_cell"))
         

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=10, #big
                             no_plot=F ,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )


plot_cnv(infercnv_obj, out_dir = "./",output_filename = paste0("infercnv",'inte7'), output_format = "pdf",color_safe_pal=F)





#================================================================================================================================
# plot infercnv score 
#================================================================================================================================

#
library(infercnv)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/B_cell_maliganat/run.final.infercnv_obj")
obs <- read.table("/public/workspace/lily/Lung2Brain/inte7/infercnv/B_cell_maliganat/infercnv.observations.txt",sep=" ")
ref <- read.table("/public/workspace/lily/Lung2Brain/inte7/infercnv/B_cell_maliganat/infercnv.references.txt",sep=" ")

#====================
as.data.frame(obs) -> obs1
obs1$chr <- dat@gene_order$chr

as.data.frame(ref) -> ref1
ref1$chr <- dat@gene_order$chr

#====================
# obs adj_normal
aggregate(.~chr,data=obs1,FUN=median) -> res_obs
rownames(res_obs) <- res_obs$chr
apply(res_obs[,-1],1,function(x){median(x)}) -> res_obs.f


#====================
# tissue normal 
aggregate(.~chr,data=ref1,FUN=median) -> res_ref
rownames(res_ref) <- res_ref$chr
apply(res_ref[,-1],1,function(x){median(x)}) -> res_ref.f

library(ggplot2)

#=========================================================================
plotdat <- as.data.frame(cbind(res_ref.f,res_obs.f))
plotdat$chr <- factor(rownames(plotdat),levels=rownames(plotdat))
tmp1 <- plotdat[,c(1,3)]
tmp2 <- plotdat[,c(2,3)]

colnames(tmp1) <- c("infercnv_score","chr")
tmp1$sample <- rep("reference",length=nrow(tmp1))
rownames(tmp1) <- NULL

colnames(tmp2) <- c("infercnv_score","chr")
tmp2$sample <- rep("malignant",length=nrow(tmp2))
rownames(tmp2) <- NULL

res <- rbind(tmp1,tmp2)


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/genome_instability.pdf",useDingbats=F)
ggplot(data = res,aes(x=chr,y=infercnv_score,color=sample,group=sample))+
	geom_line()+
	geom_point()+
	ylim(0.90,1.1)+
	theme_classic()+
	#theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
dev.off()





#================================================================================================================================================
#
# infer cnv for every sample 
#
#================================================================================================================================================
library(infercnv)
#================================================================================================================================================

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")
sample <- sample <- unique(dat$orig.ident)
# celltype <- unique(dat$type)
# celltype <- celltype[-grep("unknow",celltype)]

for(i in 1:length(sample)){
	tmp <- subset(dat,cells=which(dat$orig.ident==sample[i]&dat$type%in%c("T_cell","maliganant")))
	dir.create(paste0("/public/workspace/lily/Lung2Brain/inte7/infercnv/",sample[i]))
	setwd(paste0("/public/workspace/lily/Lung2Brain/inte7/infercnv/",sample[i]))

	write.table(tmp$type,paste(sample[i],"cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
	count<-as.matrix(tmp@assays$RNA@counts)
	write.table(count,paste(sample[i],"count_exp.txt",sep='_'),sep="\t",quote=F)

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(sample[i],"count_exp.txt",sep='_'),
         annotations_file=paste(sample[i],"cell_info.txt",sep='_'),
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names=c("T_cell"))
         

	infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=10, #big
                             no_plot=F ,
                             denoise=TRUE,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )


	plot_cnv(infercnv_obj, out_dir = "./",output_filename = paste0("infercnv",sample[i]), output_format = "pdf",color_safe_pal=F)


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=tempfile(), 
                             cluster_by_groups=TRUE, 
                             ,
                             HMM=TRUE)

}

























#===============================================================================================================================================
# CCLE 
#
#===============================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/CCLE_RNAseq_genes_rpkm_20180929.gct",skip=2,header=T)
info <- read.table("/public/workspace/lily/metastasis/data/verify/CCLE/KRAS_Cellline.txt",sep="\t",header=T)
dat.f <- dat[,c(1,2,which(colnames(dat)%in%info$CCLE.name))]
dat.f$Name <- NULL
res <- aggregate(.~Description,data=dat.f,FUN=mean)
rownames(res) <- res$Description
res$Description <- NULL

#========================= ssgsea 
#
#===============================================================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(res),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)
rownames(info) <- info$CCLE.name
res.mod <- merge(mod,info,by="row.names")

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CCLE_LUAD_KRAS.pdf")
boxplot(BMS_test_norm~KRAS_Mut,data=res.mod)
dev.off()












#=======================================================================================================================================================================================
# TCGA Methylation 
#
#=======================================================================================================================================================================================

load('~/metastasis/data/verify/TCGA_LUAD/Methylation/Methylation.RData')
load('/public/workspace/lily/Lung2Brain/inte7/TCGA_LUAD_mod.RData')
rownames(mey) <-mey[,1]
#
luad_mod[which(rownames(luad_mod)%in%colnames(mey)),] -> mod
mey[,which(colnames(mey)%in%rownames(luad_mod))] -> mey
mod <- mod[order(rownames(mod)),]
mey <- mey[,order(colnames(mey))]
#
lows <- unname(which(mod[,2]<unname(quantile(mod[,2],0.33))))
highs <- unname(which(mod[,2]>unname(quantile(mod[,2],0.67))))
#caculate 
res <- apply(mey,1,function(x){
	x<- as.numeric(x)
	c(mean(x[highs]),mean(x[lows]),(mean(x[highs])/mean(x[lows])),wilcox.test(x[highs],x[lows])$p.value)
	})

res <- t(res)
colnames(res) <- c('mean_H','mean_L','FC','p')
as.data.frame(res) -> res
res$p.adj <- p.adjust(res$p)

#=========================================
# annoation data
#
#=========================================
load('~/metastasis/data/verify/TCGA_LUAD/Methylation/ann.RData')
as.data.frame(ann) -> ann
res$Name <- rownames(res)
merge(res,ann,by='Name')-> f_res
save(f_res,file='/public/workspace/lily/Lung2Brain/inte7/res_Data/LUAD_Methy_diff.RData')

#=========================================
tmp1 <- f_res[which(f_res$p.adj<0.01&f_res$FC<0.5),]
tmp2 <- f_res[which(f_res$p.adj<0.01&f_res$FC>2),]
#a<- a[-which(a$UCSC_RefGene_Name==''),]
tmp -> luaddiffmey
save(luaddiffmey,file='./LUAD_diff_mey.RData')











#===========================================================================================================================================================================================
# map to gene mutation 
#
#===========================================================================================================================================================================================
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/vaf_luad.RData")
tmp <- as.data.frame(luad_mod)
tmp$Sample_ID <- rownames(tmp)
vafdat <- merge(vaf_luad,tmp,by="Sample_ID")
vafdat[,grep("inte11",colnames(vafdat))] <- NULL # delete old version BMS score
#============================ gene order 


aggregate(BMS_test_norm~gene,data=vafdat,FUN=mean) -> gene.detail
gene.detail <- gene.detail[order(gene.detail$BMS_test_norm),]
sample.detail <- data.frame(table(vafdat$gene))
colnames(sample.detail) <- c("gene","sample_num")
res.vaf <- merge(gene.detail,sample.detail,by='gene')
res.vaf <- res.vaf[order(res.vaf$BMS_test_norm),]

#=================================
# check if gene is mutual exclusive 
vafdat$young_type <- 0
sample1 <- unique(vafdat[which(vafdat$gene=="KEAP1"),"Sample_ID"])
vafdat$young_type[which(vafdat$Sample_ID%in%sample1)] <- 1

vafdat$old_type <- 0
sample2 <- unique(vafdat[which(vafdat$gene=="MET"),"Sample_ID"])
vafdat$old_type[which(vafdat$Sample_ID%in%sample2)] <- 1

tmp <- unique(vafdat[which(vafdat$young_type+vafdat$old_type>=1),c("Sample_ID","young_type","old_type")])
dim(tmp)
which(tmp[,2]+tmp[,3]==2)
length(which(tmp[,2]+tmp[,3]==2))/nrow(tmp)

#=========================== filter some genes 
#
#=================================================================================
#
res.vaf.f <- res.vaf[which(res.vaf$sample_num>20),]
load("/public/workspace/lily/MSK_LUAD/TMB/gene.rda")
res.vaf.ff <- res.vaf.f[which(res.vaf.f$gene%in%a),]


# famous genes 
fgene <- c("TP53","KRAS",'KEAP1','EGFR','STK11','NF1','SMARCA4','BRAF','CDKN2A','MET','RBM10','PIK3CA','RIT1','U2AF1',
	'ATM','ARID1A','SLC4A5','RB1','ERBB2','MAP2K1','STX2','NBPF1',
	'MLL3','FAT1','APC','ARID2','SLC39A6','ARHGP35','SMAD4','CTNNB1','CDK12','MBD1','FBXO16','PBLD','ZNF774','NRAS','CDKN1B')

res.vaf.ff <- res.vaf[which(res.vaf$gene%in%fgene),]


mskdat <- read.delim("~/MSK_LUAD/TMB/data_mutations_mskcc.txt",sep="\t",header=T)
mskdat$MAF <- mskdat$t_alt_count/(mskdat$t_ref_count+mskdat$t_alt_count) # caculaye maf 
mskdat.f <- mskdat[,c(colnames(mskdat)[1:20],"MAF")] #filter some information

#======================================================
youngsample <- unique(mskdat.f[which(mskdat.f$Hugo_Symbol%in%c("STK11","KEAP1","RB1")),"Tumor_Sample_Barcode"])
oldsample <- unique(mskdat.f[which(mskdat.f$Hugo_Symbol%in%c("MET","ERBB2","FAT1")),"Tumor_Sample_Barcode"])

#======================================================
# some sample have old and young mutation 

tmp <- intersect(youngsample,oldsample)
if(length(tmp)>0){
	# for(i in 1:length(tmp){
	# 	tmp_sample <- mskdat.f[which(mskdat.f$Tumor_Sample_Barcode==tmp[1]),]
	# }

	#=========================================================
	# do not know how to handle just delete 
	youngsample <- youngsample[-which(youngsample%in%tmp)]
	oldsample <- oldsample[-which(oldsample%in%tmp)]

}


#===============================================================================================================================================================================clinical dat
mskclin <- read.delim("~/MSK_LUAD/TMB/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin[which(mskclin$CANCER_TYPE_DETAILED=="Lung Adenocarcinoma"),c(1,2)] -> luad_sample # all LUAD patients
#==================================================================================
ysample.f <- youngsample[which(youngsample%in%luad_sample$SAMPLE_ID)]
osample.f <- oldsample[which(oldsample%in%luad_sample$SAMPLE_ID)]
both <- oldsample[which(tmp%in%luad_sample$SAMPLE_ID)]

# #=============================================================
# mskdat.f.luad <- mskdat.f[which(mskdat.f$Tumor_Sample_Barcode%in%luad_sample$SAMPLE_ID),] # mutation info for luad patient
# save(mskdat.f.luad,file="./MSK_LUAD/TMB/mskdat.f.luad.mut.RData")

#====================survival dat
msksurv <- read.delim("~/MSK_LUAD/TMB/data_clinical_patient.txt",sep="\t",comment.char="#")
substr(ysample.f,1,9) ->ysample.f
substr(osample.f,1,9) ->osample.f  #change names to map 
substr(both,1,9) -> both.f
msksurv.f <- msksurv[which(msksurv$PATIENT_ID%in%c(osample.f,ysample.f,both.f)),]

#====================
# run surv analysis 
library(survminer)
library(survival)

msksurv.f$type <- "BMS_low"
msksurv.f$type[which(msksurv.f$PATIENT_ID%in%osample.f)] <- "BMS_high"
msksurv.f$type[which(msksurv.f$PATIENT_ID%in%both.f)] <- "BMS_low"
msksurv.f$OS_status <- 0 #change type info 
msksurv.f$OS_status[which(msksurv.f$OS_STATUS=="DECEASED")] <- 1

table(msksurv.f$type)
surv <- Surv(msksurv.f$OS_MONTHS,msksurv.f$OS_status)
km <- survfit(surv~type,data=msksurv.f)
ggsurvplot(km,pval=T)










#========================================================================================================================================================
# plot point 
#
#=========================================================================================================================
res.vaf.ff$type <- "low_BMS"
res.vaf.ff$type[which(res.vaf.ff$BMS_test_norm>median(res.vaf.ff$BMS_test_norm))] <- "high_BMS"
res.vaf.ff$label=""
res.vaf.ff$label[which(res.vaf.ff$sample_num>=20)] <- as.vector(res.vaf.ff$gene[which(res.vaf.ff$sample_num>=20)])
ggplot(res.vaf.ff, aes(x=BMS_test_norm, y=sample_num, size=sample_num,color=type))+ geom_point()+
	geom_hline(aes(yintercept=20))+
	geom_vline(aes(xintercept=median(res.vaf.ff$BMS_test_norm)))+
	geom_text_repel(aes(label=res.vaf.ff$label),check_overlap = TRUE, size=3.5)+theme_classic()






#========================================================================================================================================================
# 
# GSE68465
library(idmap1,lib.loc="/public/workspace/wangzy/R/x86_64-pc-linux-gnu-library/3.6.0")
ids=getIDs('gpl96')
mat <- read.table('/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_series_matrix.txt',sep="\t",comment.char="!",header=T)
info <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS")
as.data.frame(info) -> info

colnames(mat)[1] <- "probe_id"
mat.merge <- merge(ids,mat,by="probe_id")
mat.merge$probe_id <- NULL
mat.merge$gpl <- NULL
aggregate(.~symbol,data=mat.merge,FUN=mean) -> res

rownames(res) <- res$symbol
res$symbol <- NULL


#================================================
# run ssgsea to caculate 
#================================================
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
pdf("./GSE68465.pdf")
ggsurvplot(km,pval=T,surv.median.line='hv')$plot+scale_colour_manual(values=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'))
dev.off()









#==================================================================================================================================
# GSE14995 analysis JAG-1 invasion 
#
#==================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14995/GSE14995_series_matrix.txt",sep="\t",header=T,comment.char="!")
library(idmap1,lib.loc="/public/workspace/wangzy/R/x86_64-pc-linux-gnu-library/3.6.0")
ids=getIDs('gpl570')
GPL570 <- ids[,c("probe_id","symbol")]
colnames(GPL570)[1] <- "ID_REF"
dat.f <- merge(GPL570,dat,by="ID_REF")
#======================== prepare 
dat.f$ID_REF <- NULL
aggregate(.~symbol,data=dat.f,FUN=mean) -> res
rownames(res) <- res$symbol
res$symbol <- NULL
saveRDS(res,file="/public/workspace/lily/metastasis/data/verify/GSE14995/GSE14995_dat.RDS")
#======================== run ssgsea 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(res),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)
mod <- as.data.frame(mod)
mod$type <- "low_invasion"
mod$type[4:6] <- "high_invasion"
#======================== plot result 
boxplot(BMS_test_norm~type,data=mod)





#=======================================================================================================================================
# GSE123066
# gefitinib ressitance 
#
#=======================================================================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE123066/GSE123066_series_matrix.txt",sep="\t",header=T,comment.char="!")
rownames(dat) <- dat$ID_REF
dat$ID_REF <- NULL
#========================== run ssgsea 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)
mod <- as.data.frame(mod)





















#=================================================================================================================================================================
# integrated GBM data 
#==============================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/inte7_ann.RDS")
GBM1 <- readRDS("/public/workspace/lily/PS/Final_716/lesion1.RDS")
GBM2 <- readRDS("/public/workspace/lily/PS/Final_716/lesion2.RDS")

integration.anchors <- FindIntegrationAnchors(object.list = c(dat,GBM1,GBM2))
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

#=================================================================================================================================
# add some information 
#
#=================================================================================================================================
inte$type_group[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- "GBM" # sample info 

# cell type info 
inte$type[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))] <- inte$celltype[which(inte$orig.ident%in%c('RD-20180817-001-SR18271','RD-20180817-002-SR18271'))]
inte$type[which(inte$type%in%c("BMDM","MG"))] <- "Myeloid"
inte$type[which(inte$type=="Tumor_cell")] <- "maliganant"

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")







#===============================================================================================================================
#
# Treg 
#================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/Tcell_recluster.f.RDS")
DefaultAssay(dat) <- "RNA"
FeaturePlot(dat,features=c(""))


















#=============================================================================================================================
# check TCGA surv and MSK surv 
# 
#=============================================================================================================================
# TCGA mut 
load("/public/workspace/lily/Lung2Brain/inte7/TCGA_LUAD_mod.RData")
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/vaf_luad.RData")
tmp <- as.data.frame(luad_mod)
tmp$Sample_ID <- rownames(tmp)
vafdat <- merge(vaf_luad,tmp,by="Sample_ID")
vafdat[,grep("inte11",colnames(vafdat))] <- NULL # delete old version BMS score
#========================= get sample 
unique(vafdat[which(vafdat$gene%in%c("MET","ERBB2","FAT1")),"Sample_ID"]) -> tcga.sample
# get survival time 
#=========================
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
ann -> luad_ann
rm(ann)
luad_ann[which(gsub("-",".",rownames(luad_ann))%in%tcga.sample),] -> surv.tcga 
surv.tcga$OS.month <- surv.tcga$OS.time/30 # change to months

#=================
# MSK 

msksurv.f[which(msksurv.f$type=="BMS_high"),] -> surv.msk

#=================
# combine
#=================
surv.tcga[,c("OS","OS.month")] -> surv.tcga.dat
surv.tcga.dat$group <- "TCGA"
colnames(surv.tcga.dat)[1:2] <- c("OS","OS.month")
surv.msk[,c("OS_status","OS_MONTHS")] -> surv.msk.dat
surv.msk.dat$group <- "MSK"
colnames(surv.msk.dat)[1:2] <- c("OS","OS.month")

res.f <- rbind(surv.tcga.dat,surv.msk.dat)


library(survminer)
library(survival)

surv <- Surv(res.f$OS.month,res.f$OS)
km <- survfit(surv~group,data=res.f)

ggsurvplot(km,pval=T)



#======================================================================================================
# check if TCGA sample and MSK sample have difference in survival 
#
#======================================================================================================

load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
ann -> luad_ann
rm(ann)
luad_ann[,c("OS","OS.time")] -> surv.tcga 
surv.tcga$OS.month <- surv.tcga$OS.time/30 # change to months
surv.tcga$OS.time=NULL
surv.tcga$group="TCGA"

mskclin <- read.delim("~/MSK_LUAD/TMB/data_clinical_sample.txt",sep="\t",header=T,comment.char="#")
mskclin[which(mskclin$CANCER_TYPE_DETAILED=="Lung Adenocarcinoma"),c(1,2)] -> luad_sample # all LUAD patients
msksurv <- read.delim("~/MSK_LUAD/TMB/data_clinical_patient.txt",sep="\t",comment.char="#")
msksurv.f <- msksurv[which(msksurv$PATIENT_ID%in%luad_sample$PATIENT_ID),]

msksurv.f[,c("OS_STATUS","OS_MONTHS")] -> surv.msk
as.numeric(surv.msk$OS_STATUS) -> surv.msk$OS_STATUS
surv.msk$OS_STATUS[which(surv.msk$OS_STATUS==2)] <- 0
surv.msk$group="MSK"
colnames(surv.msk)[1:2] <- c("OS","OS.month")

res.f <- rbind(surv.tcga,surv.msk)



#=================================================================================================================================================================
# check TCGA mutation sample 
#=======================================================
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
ann -> luad_ann
rm(ann)
luad_ann[,c("OS","OS.time")] -> surv.tcga 
rownames(surv.tcga) <- gsub("-",".",rownames(surv.tcga))

load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/vaf_luad.RData")
dat.high <- as.vector(unique(vaf_luad[which(vaf_luad$gene%in%c("SPTBN1","SLC9A2","GALNT13")),"Sample_ID"]))
dat.low <- as.vector(unique(vaf_luad[which(vaf_luad$gene%in%c("STK11","KEAP1","RB1")),"Sample_ID"]))

tmp <- intersect(dat.high,dat.low)
dat.high.f <- dat.high[-which(dat.high%in%tmp)]
dat.low.f <- dat.low[-which(dat.low%in%tmp)]

 
info.surv <- surv.tcga[which(rownames(surv.tcga)%in%c(dat.high.f,dat.low.f)),]
info.surv$type <- "low"
info.surv$type[which(rownames(info.surv)%in%dat.high.f)] <- "high"


library(survminer)
library(survival)

surv <- Surv(info.surv$OS.time,info.surv$OS)
km <- survfit(surv~type,data=info.surv)

ggsurvplot(km,pval=T)



#=====================================================================================
# another ways to calculate 
#
#=====================================================================================
load('/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_anno.RData')
ann -> luad_ann
rm(ann)
luad_ann[,c("OS","OS.time")] -> surv.tcga 
rownames(surv.tcga) <- gsub("-",".",rownames(surv.tcga))

load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/vaf_luad.RData")
res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(vaf_luad,res.vaf,by="gene")
res.f <- aggregate(BMS_test_norm~Sample_ID,data=res,FUN=mean)
#========================================================================================================
surv.tcga$type <- "Median"
surv.tcga$type[which(rownames(surv.tcga)%in%res.f$Sample_ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.67))])] <- "High"
surv.tcga$type[which(rownames(surv.tcga)%in%res.f$Sample_ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.33))])] <- "Low"

library(survminer)
library(survival)

surv <- Surv(surv.tcga$OS.time,surv.tcga$OS)
km <- survfit(surv~type,data=surv.tcga)

ggsurvplot(km,pval=T)





#=================================================================================================================================================================
# MSK35 sample 
#
#=================================================================================================================================================================
dat <- read.table("/public/workspace/lily/MSK_LUAD/MSK_35/msk35_mutation.txt",sep="\t",header=T) # mutation data 
info <- read.table("/public/workspace/lily/MSK_LUAD/MSK_35/msk35_sampleinfo.txt",sep="\t",header=T)



res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(dat,res.vaf,by="gene")

res.f <- aggregate(BMS_test_norm~Sample,data=res,FUN=median)

dat.high <- res.f$Sample[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$Sample[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]





#===============================================================
# 
dat.high <- as.vector(unique(dat[which(dat$Gene%in%c("MET","INSRR","PDGFRA")),"Sample"]))
dat.low <- as.vector(unique(dat[which(dat$Gene%in%c("STK11","KEAP1","APC")),"Sample"]))

#===============================================================
# 
tmp <- c(as.vector(dat.high),as.vector(dat.low))
info.surv <- info[which(as.vector(info$Study_ID)%in%tmp),c("Study_ID","PFS","Event","Resp","Durable_Clinical_Benefit")]
info.surv$type <- "low"
info.surv$type[which(info.surv$Study_ID%in%dat.high)] <- "high"


#===============================================================
# surv plot
#===============================================================
library(survminer)
library(survival)

surv <- Surv(info.surv$PFS,info.surv$Event)
km <- survfit(surv~type,data=info.surv)





#====================================================================================================================================================================
# GSE126044
# counts to TPM
#
#====================================================================================================================================================================
count <- read.table("/public/workspace/lily/metastasis/data/verify/GSE126044/GSE126044_counts.txt",sep="\t",header=T)
colnames(count)[1] <- "Gene"
gene.info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",sep="\t",header=T)
colnames(gene.info) <- c("Gene","Chr","start","end")
gene.info$Len <- gene.info$end-gene.info$start
#========================================================================
# just 18253 genes can find Length 
#
#========================================================================

count.f <- merge(count,gene.info[,c("Gene","Len")],by="Gene")
apply(count.f[,-1],1,function(x){x[]/x[length(x)]}) -> tmp
colnames(tmp)<- count.f$Gene
res <- apply(tmp,1,function(x){x[]/sum(x)})
res.f <- res*(10^6)
res.ff <- res.f[,-grep("Len",colnames(res.f))]

saveRDS(res.ff,file="/public/workspace/lily/metastasis/data/verify/GSE126044/TPM.RDS")



#=========================================================================
# calculate ssgsea
# 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE126044/TPM.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),"BMS_test","/public/workspace/lily/Lung2Brain/inte7/",permN=1000)

mod <- data.frame(mod)
mod$type <- "non-responder"
mod$type[which(rownames(mod)%in%c("Dis_02","Dis_15","Dis_04","Dis_17","Dis_10"))] <- "responder"
aggregate(BMS_test_norm~type,data=mod,FUN=median)
















#==========================================================================================================
# GSE135222
# prepare 
#
#==========================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/surv_info.txt",sep="\t")
dat <- t(dat)
colnames(dat) <- c("Status","PFS","ID")
dat.f <- dat[-1,]
sapply(strsplit(dat.f[,1],": "),function(x){x[2]}) -> dat.f[,1]
sapply(strsplit(dat.f[,2],": "),function(x){x[2]}) -> dat.f[,2]
#====================================
rownames(dat.f)<- dat.f[,"ID"]
dat.f <- data.frame(dat.f)
as.numeric(as.vector(dat.f$Status)) -> dat.f$Status
as.numeric(as.vector(dat.f$PFS)) -> dat.f$PFS
saveRDS(dat.f,file="/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.RDS")

#==============================================================================================
# transform gene 
#
#==============================================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(Seurat)

tmp <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/GSE135222_GEO_RNA-seq_omicslab_exp.tsv",sep="\t",header=T)
sapply(strsplit(as.vector(tmp$gene_id),"\\."),function(x){x[1]}) -> tmp$ENSEMBL
gene_name <- bitr(as.vector(tmp$ENSEMBL), fromType="ENSEMBL", toType=c("SYMBOL", "GENENAME"), OrgDb="org.Hs.eg.db")

res.f <- merge(tmp,gene_name[,1:2],by="ENSEMBL")
res.f$ENSEMBL <- NULL
res.f$gene_id <- NULL

res.ff <- aggregate(.~SYMBOL,data=res.f,FUN=mean)
rownames(res.ff) <- res.ff$SYMBOL
res.ff$SYMBOL <- NULL
saveRDS(res.ff,file="/public/workspace/lily/metastasis/data/verify/GSE135222/expr.RDS")


#===================================================================================================================================================================================================
# calculate ssgsea
#
#=================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),"BMS_test","/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- data.frame(mod)


#=================================================================================================================================
# info 
sample.info <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.RDS")
sample.info$type <- "Median"
sample.info$type[which(mod$BMS_test_norm<median(mod$BMS_test_norm))] <- "Low"
sample.info$type[which(mod$BMS_test_norm>median(mod$BMS_test_norm))] <- "high"


library(survminer)
library(survival)

surv <- Surv(sample.info$PFS,sample.info$Status)
km <- survfit(surv~type,data=sample.info)

ggsurvplot(km,pval=T)








#=================================================================================================================================================================================
# Cancer cell immunotherapy
#
#=================================================================================================================================================================================

dat <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/mutation.txt",sep="\t",header=T)
sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/cancercell/sample.txt",header=T,sep="\t")
sample.info.f <- sample.info[which(sample.info$Histology=="non-squamous"),] # non-squamous 

# 
dat.high <- as.vector(unique(dat[which(dat$gene_name%in%c("SPTBN1","SLC9A2","GALNT13")),"Patient_ID"]))
dat.low <- as.vector(unique(dat[which(dat$gene_name%in%c("ZNF592","APBA2","ATR")),"Patient_ID"]))
tmp <- intersect(dat.high,dat.low)
dat.high.f <- dat.high  [-which(dat.high%in%tmp)]
dat.low.f <- dat.low   [-which(dat.low%in%tmp)]




# 
info.surv <- sample.info.f[which(sample.info.f$Patient.ID%in%c(dat.high.f,dat.low.f)),c("Patient.ID","PFS..months.","PFS...0.censor..1.event.","Best.Overall.Response","Clinical.benefit..DCB...durable.clinical.benefit..NDB...no.durable.benefit")]
info.surv$type <- "low"
info.surv$type[which(info.surv$Patient.ID%in%dat.high.f)] <- "high"
colnames(info.surv) <- c("Patient.ID","PFS.time","PFS","Response","DCB","type")


library(survminer)
library(survival)

surv <- Surv(info.surv$PFS.time,info.surv$PFS)
km <- survfit(surv~type,data=info.surv)

ggsurvplot(km,pval=T)







#=============================================================================================================================================

dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)

ann <- dat[,c(1,31,32,44,45,52,53,125)]
ann$Sample_ID <- gsub("-",".",ann$sampleID)

load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/vaf_luad.RData")
dat.high <- as.vector(unique(vaf_luad[which(vaf_luad$gene%in%c("MET","ERBB2","FAT1")),"Sample_ID"]))
dat.low <- as.vector(unique(vaf_luad[which(vaf_luad$gene%in%c("STK11","KEAP1","RB1")),"Sample_ID"]))


tmp <- intersect(dat.high,dat.low)
if(length(tmp)>0){
	dat.high.f <- dat.high[-which(dat.high%in%tmp)]
	dat.low.f <- dat.low[-which(dat.low%in%tmp)]

}else{
	dat.low.f <- dat.low
	dat.high.f <- dat.high
}


#===========================================================================================================================================
# 
ann.f <- ann[-which(ann$targeted_molecular_therapy)]

surv.ann <- ann.f[which(ann.f$Sample_ID%in%c(dat.low.f,dat.high.f)),]
surv.ann$type <- "unknow"
surv.ann$type[which(surv.ann$Sample_ID%in%dat.high.f)] <- "high"
surv.ann$type[which(surv.ann$Sample_ID%in%dat.low.f)] <- "low"


library(survminer)
library(survival)

surv <- Surv(surv.ann$OS.time,surv.ann$OS)
km <- survfit(surv~type,data=surv.ann)





#=======================================================================================================================================
# try to classifiy LUAD and LUSC samples 
#===================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/60sample_mutation.txt",sep="\t",header=T)
sapply(strsplit(gsub("\\(",":",as.vector(dat$Gene)),":"),function(x){x[1]}) -> dat$gene

dat.high <- as.vector(unique(dat[which(dat$gene%in%c("MET","ERBB2","FAT1")),"Sample.ID"]))
dat.low <- as.vector(unique(dat[which(dat$gene%in%c("STK11","KEAP1","RB1")),"Sample.ID"]))

#===================================================
# surv info 
#===================================================
sample.info <- read.table("/public/workspace/lily/metastasis/data/verify/GSE135222/sample_info.txt",sep="\t",header=T)

sample.info[which(sample.info$Patient.ID%in%c(dat.high,dat.low)),] -> surv.info
surv.info$type <- "low"
surv.info$type[which(surv.info$Patient.ID%in%c(dat.high))] <- "high"

#=============================

library(survminer)
library(survival)

surv <- Surv(surv.info$PFS,surv.info$PFS_event)
km <- survfit(surv~type,data=surv.info)






res.vaf <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res.vaf.RDS")
res <- merge(dat,res.vaf,by="gene")

res.f <- aggregate(BMS_test_norm~Sample.ID,data=res,FUN=median)

dat.high <- res.f$Patient_ID[which(res.f$BMS_test_norm>quantile(res.f$BMS_test_norm,0.5))]
dat.low <- res.f$Patient_ID[which(res.f$BMS_test_norm<quantile(res.f$BMS_test_norm,0.5))]









#============================================================================================
# TCGA purity
#============================================================================================
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/LUAD_purity.RData")
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")

mod <- data.frame(TCGAlong.id=rownames(luad_mod),score=luad_mod[,2])
pur <- data.frame(TCGAlong.id=gsub('-','.',substr(luad_purity$tcga_id,1,15)),pur=luad_purity$abs_purity)


res <- merge(mod,pur,by='TCGAlong.id')
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_LUAD_purity.pdf",useDingbats=F)
ggplot(data =res,aes(x =score,y =pur)) + 
  geom_point(colour = "#426671", size = 2) +  
  geom_smooth(method = lm,colour='#764C29',fill='#E7E1D7')+
  stat_cor(method = "pearson")+
  theme_classic()+
  labs(x="BMS",y="Purity",title='LUAD')
dev.off()





#=======================================================================================
# GBM plot
#
#=======================================================================================
library(Seurat)
lesion1 <- readRDS("/public/workspace/lily/PS/data/hbrd001.rds")
lesion2 <- readRDS("/public/workspace/lily/PS/data/hbrd002.rds")


pdf("./GBM_landscape.pdf",useDingbats=F)
DimPlot(tmp,group.by="llymarker",cols=c('#910241',"#984EA3","#F29403","#8CA77B","#B2DF8A","#E41A1C"))
barplot(as.matrix(as.numeric(as.vector(table(tmp$llymarker)))/21518),beside=F,col=c('#910241',"#984EA3","#F29403","#8CA77B","#B2DF8A","#E41A1C"))
dev.off()


#          B_cell     Endothelial      Fibroblast      maliganant         Myeloid
#            1884             392            1046           10321            9429
# Oligodendrocyte          T_cell          unknow
#             954           11874             180
cols <- c('#377EB8','#910241','#984EA3','#E41A1C','#F29403','#8CA77B','','#999999')

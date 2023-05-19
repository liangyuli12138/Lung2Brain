
#===========================================================================================
# 2020-10-5
# Lung2Brain Figure1 
#===========================================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 

#          B_cell     Endothelial      Fibroblast      maliganant         Myeloid
#            1884             392            1046           10321            9429
# Oligodendrocyte          T_cell          unknow
#             954           11874             180
cols <- c('#377EB8','#910241','#984EA3','#E41A1C','#F29403','#8CA77B','#B2DF8A','#999999')

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/lanscape.pdf",useDingbats=F)
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
# cell type percentage check 
#==============================================================================================================
table(dat$type,dat$type_group) -> tmp
tmp.f <- tmp[-which(rownames(tmp)=="maliganant"),]
apply(tmp.f,2,function(x){x/sum(x)})-> tmp.res
# maybe should filter maligant cells 
cols <- c('#377EB8','#910241','#984EA3','#F29403','#8CA77B','#B2DF8A','#999999')
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/celltype_percentage.pdf")
barplot(tmp.res,col=cols)
legend("topright",title="Cell component", c("B cell","Endothelial","Fibroblast","Myeloid","Oligodendrocyte","T cell","unknow"),pch=c(15,15,15,15,15,15,15), col=cols)
dev.off()














#==============================================================================================================
# expression genes between tumor/non-tumor and LC/LCBM
# genes 
#==============================================================================================================
data <- dat[['RNA']]@data

tumor_gene <- apply(data[,which(dat$maliganant=="tumor")],2,function(x){length(which(x[]>0))})
ntumor_gene <- apply(data[,which(dat$maliganant=="non-tumor")],2,function(x){length(which(x[]>0))})

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/tumor_ntumor_gene.pdf")
boxplot(tumor_gene,ntumor_gene,names=c('malignant','non-malignant'),main="express gene numbers",col=c("#E41A1C","#A6CEE3"),outline=F)
legend('topright',legend=paste('p =',round(wilcox.test(tumor_gene,ntumor_gene)$p.value,5)),bty='n')
dev.off()


#=====================================================
# compare LUAD and LCBM tumor cell gene expression numbers
#=====================================================
data <- dat[['RNA']]@data
LC_gene <- apply(data[,which(dat$maliganant=="tumor"&dat$type_group=="LC")],2,function(x){length(which(x[]>0))})
LCBM_gene <- apply(data[,which(dat$maliganant=="tumor"&dat$type_group=="LCBM")],2,function(x){length(which(x[]>0))})
# the differencr maybe caused by batch effect 
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/tumorgene_LC_LCBM.pdf")
boxplot(LC_gene,LCBM_gene,names=c('LC_tumor','LCBM_tumor'),main="express gene numbers",col=c("#E41A1C","#A6CEE3"),outline=F)
legend('topright',legend=paste('p =',round(wilcox.test(LC_gene,LCBM_gene)$p.value,5)),bty='n')
dev.off()






#===============================================================================================================
# stemness calculate 
#
#===============================================================================================================
conda activate CytoTrace
R
# cytotrace should use R-4.0.2
#==========================
library(Seurat)
library(CytoTRACE)


dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
tumor <- subset(dat,cells=which(dat$maliganant=="tumor"))
results <- CytoTRACE(as.matrix(tumor[["RNA"]]@data),ncores = 4,subsamplesize = 1000)

tumor$cytotrace <- results$CytoTRACE
save(results,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/cytotrace_result.RData")

#=========================
# plot result 
#=========================
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/tumor_cytotrace.pdf")
boxplot(cytotrace~type_group,data=tumor@meta.data,outline=F)
legend('topright',legend=paste('p =',round(wilcox.test(cytotrace~type_group,data=tumor@meta.data)$p.value,5)),bty='n')
dev.off()

library(ggplot2)
ggplot(tumor@meta.data, aes(x = type_group,y=cytotrace))+ geom_violin()

# try to plot DimPlot 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
tumor <- subset(dat,cells=which(dat$maliganant=="tumor"))
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/cytotrace_result.RData")
tumor$cytotrace <- results$CytoTRACE
FeaturePlot(tumor,features="cytotrace")

#=====================================
# try to find same result in Samgsung
#=====================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")
dat <- CellCycleScoring(dat ,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
sub.tmp <- subset(dat,cells=which(dat$Cell_type.refined=="Epithelial cells"))
sub.dat <- subset(sub.tmp,cells=which(sub.tmp$type%in%c("Brian_metastasis","tumor_Lung_advanced","tumor_Lung_early")))

#================================================
# use cell cycle do not have good ways to reslove 
# try to use cytotrace 
#================================================
library(Seurat)
library(CytoTRACE)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")
sub.tmp <- subset(dat,cells=which(dat$Cell_type.refined=="Epithelial cells"))
sub.dat <- subset(sub.tmp,cells=which(sub.tmp$type%in%c("Brian_metastasis","tumor_Lung_advanced","tumor_Lung_early")))
results <- CytoTRACE(as.matrix(sub.dat[["RNA"]]@data),ncores = 4,subsamplesize = 1000)
sub.dat$cytotrace <- results$CytoTRACE

#================================================
# 












#=====================================================================================================================
# ROGUE calculate 
# purity 
#=====================================================================================================================
library(rlist)
library(Seurat)
library(pheatmap)
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))


dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
tumor <- subset(dat,cells=which(dat$maliganant=="tumor"))

LC_tumor <- subset(tumor,cells=which(tumor$type_group=="LC"))
LCBM_tumor <- subset(tumor,cells=which(tumor$type_group=="LCBM"))

#================================================================
# Lung Cancer
#################################################################
expr_LC <- as.matrix(LC_tumor[["RNA"]]@data)
expr_LC <- matr.filter(expr_LC, min.cells = 10, min.genes = 10) #filter gene and cells 
ent.res.LC <- SE_fun(expr_LC)
rogue.value.LC <- CalculateRogue(ent.res.LC, platform = "UMI")


#================================================================
# Brian Metastasis
expr_LCBM <- as.matrix(LCBM_tumor[["RNA"]]@data)
expr_LCBM <- matr.filter(expr_LCBM, min.cells = 10, min.genes = 10) #filter gene and cells 
ent.res.LCBM <- SE_fun(expr_LCBM)
rogue.value.LCBM <- CalculateRogue(ent.res.LCBM, platform = "UMI")


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LC_LCBM_tumor_ROGUE.pdf")
barplot(c(rogue.value.LC,rogue.value.LCBM),main="ROGUE score",names.arg = c("LC","LCBM"))
dev.off()







#=======================================================================================================================
# GSE131907
# tmp file 
#=======================================================================================================================
dat$type <- "unsure"
dat$type[which(dat$Sample%in%c("LN_01","LN_02","LN_03","LN_04","LN_05","LN_06","LN_07","LN_08","LN_11","LN_12"))] <- "normal_Lymph_Node"

dat$type[which(dat$Sample%in%c("LUNG_N01","LUNG_N06","LUNG_N08","LUNG_N09","LUNG_N18",
				"LUNG_N19","LUNG_N20","LUNG_N28","LUNG_N30","LUNG_N31","LUNG_N34"))] <- "normal_Lung"

dat$type[which(dat$Sample%in%c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19",
				"LUNG_T20","LUNG_T25","LUNG_T28","LUNG_T30","LUNG_T31","LUNG_T34"))] <- "tumor_Lung_early"

dat$type[which(dat$Sample%in%c("NS_02","NS_03","NS_04","NS_06","NS_07","NS_12","NS_13","NS_16","NS_17","NS_19"))] <- "Brian_metastasis"

dat$type[which(dat$Sample%in%c("EFFUSION_06","EFFUSION_11","EFFUSION_12","EFFUSION_13","EFFUSION_64"))] <- "pleural_effusion"

dat$type[which(dat$Sample%in%c("EBUS_06","EBUS_28","EBUS_49","BRONCHO_58"))] <- "tumor_Lung_advanced"

dat$type[which(dat$Sample%in%c("EBUS_10","BRONCHO_11","EBUS_12","EBUS_13","EBUS_15","EBUS_19","EBUS_51"))] <- "Metastasis_Lymph_node"



LC_gene <- apply(data[,which(tumor$type=="tumor_Lung_early")],2,function(x){length(which(x[]>0))})
LCBM_gene <- apply(data[,which(tumor$type=="Brian_metastasis")],2,function(x){length(which(x[]>0))})




#=====================================================================================================================
# T cell different 
#=====================================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
tcell <- subset(dat,cells=which(dat$type=="T_cell"))

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/T_cell_fig1.pdf",useDingbats=F)
DimPlot(tcell,group.by="type_group",cols=c('#54B0E4','#F29403'),pt.size=1)
barplot(c(8030/13332,3844/22748),names=c("Lung cancer","Brian metastasis"),main="T cell percentages",col=c('#54B0E4','#F29403'))
dev.off()



#======================================================================================================================
# Meyloid cell difference
#
#=====================================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
Myeloid <- subset(dat,cells=which(dat$type=="Myeloid"))

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Myeloid_fig1.pdf",useDingbats=F)
DimPlot(Myeloid,group.by="type_group",cols=c('#54B0E4','#F29403'),pt.size=0.8)
barplot(c(1252/13332,8177/22748),names=c("Lung cancer","Brian metastasis"),main="Myeloid percentages",col=c('#54B0E4','#F29403'))
dev.off()







#====================================================================================================================
#  marker plot
#
#====================================================================================================================

marker <- c("CD3D","CD3E","CD2","PTPRC",
			"CD19","CD79A","MS4A1",
			"CD68","FCGR1A","FCGR3A",
			"MAG","MOG","CNDP1",
			"COL1A2","THY1","DCN","CD248",
			"CLDN5","CDH5","VWF","ABCG2",
			"EGFR","EPCAM","CDH1"
			)

#==========================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
data <- as.matrix(dat[["RNA"]]@data)

data.f <- data[marker,]
data.f <- t(data.f)
data.f <- data.frame(data.f)
data.f$type <- dat$type

#===========================================
# mean expression
#===========================================
aggregate(.~type,data=data.f,FUN=mean) -> mean_exp
expr.dat <- melt(mean_exp,id.vars="type")
colnames(expr.dat) <- c("cell_type","marker","mean_exp")
#==========================================
# fraction 
#==========================================
data.f.bak <- data.f
data.f.bak[data.f.bak[]>0] <- 1
data.f.bak$type <- dat$type
aggregate(.~type,data=data.f.bak,FUN=sum) -> frac
rownames(frac) <- frac[,1]
frac <- frac[,-1]
frac$cell_num <- as.numeric(table(dat$type))
# melt(frac,id="type")
res <- apply(frac,1,function(x){scale(x/x[25],center=F)})
res.f <- res[-25,]
rownames(res.f) <- colnames(frac)[1:24]
frac.dat <- melt(res.f,id="col.names")
colnames(frac.dat) <- c("marker","celltype","frac")

#=========================================
# merge result 
#=========================================
frac.dat <- frac.dat[order(frac.dat$marker),]
expr.dat$marker <- as.vector(expr.dat$marker)
expr.dat <- expr.dat[order(expr.dat$marker),]
res.final <- cbind(expr.dat,frac.dat[,"frac"])
colnames(res.final)[4] <- "frac" 

#=========================================
# set factor level 
#=========================================
res.final$cell_type <- factor(res.final$cell_type,levels=c("B_cell","T_cell","Endothelial","Fibroblast","Oligodendrocyte","Myeloid",
	"maliganant","unknow"))

res.final$marker <- factor(res.final$marker,levels=c("CD79A","MS4A1","CD19","PTPRC","CD3D","CD3E","CD2","CLDN5","VWF","ABCG2","THY1","CDH5",
	"COL1A1","COL1A2","DCN","CD248","MAG","MOG","CNDP1","CD68","FCGR3A","FCGR1A","EGFR","EPCAM","CDH1"))




pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/cell_marker_bubble.pdf",useDingbats=F)
ggplot(res.final, aes(x=marker, y=cell_type, size=frac)) + geom_point(aes(colour = mean_exp))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	scale_colour_gradient(low = "white",high = "ForestGreen")+
	theme(panel.grid.major = element_blank())
dev.off()

	
	












#=============================================
# GSE131907

mod <- data.frame(mod_lungae)
mod$type <- tumor$orig.ident
mod1 <- mod[which(mod$type=="LungTumor_early"),]
mod2 <- mod[which(mod$type=="LungTumor_advanced"),]
mod1 <- mod1[order(mod1$BMS_test_norm),]
mod2 <- mod2[order(mod2$BMS_test_norm),]

mod1$order.idx <- seq(1:nrow(mod1))/nrow(mod1)
mod2$order.idx <- seq(1:nrow(mod2))/nrow(mod2)
res <- rbind(mod1,mod2)

library(ggplot2)

pdf("./GSE131907.pdf",useDingbats=F)
ggplot(data = res,aes(x=order.idx,y=BMS_test_norm,colour=type))+
	geom_point()+theme_classic()
dev.off()

	+
	#stat_boxplot(geom = "errorbar",width=0.25)+
	stat_compare_means(comparisons=my_comparisons)+ # Add pairwise 
	labs(x="Type",y="BMS",title=modname)+
	theme_classic()+
	theme(plot.title = element_text(hjust = 0.5))








#=============================================================================================
# calculate celltype pearson 
#=============================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")

LCBM.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))
LC.dat <- subset(dat,cells=which(dat$type_group=="LC"))

#=============================================================================================
# 
LCBM.expr <- data.frame(t(as.matrix(LCBM.dat[["RNA"]]@data)))
LCBM.expr$celltype <- LCBM.dat$type
LCBM.res <- aggregate(.~celltype,data=LCBM.expr,FUN=mean)
rownames(LCBM.res) <- LCBM.res$celltype
LCBM.res$celltype <- NULL
LCBM.res.f <- t(LCBM.res)
colnames(LCBM.res.f) <- paste0("LCBM.",colnames(LCBM.res.f))
#=============================================================================================
# 
LC.expr <- data.frame(t(as.matrix(LC.dat[["RNA"]]@data)))
LC.expr$celltype <- LC.dat$type 
LC.res <- aggregate(.~celltype,data=LC.expr,FUN=mean)
rownames(LC.res) <- LC.res$celltype
LC.res$celltype  <- NULL
LC.res.f <- t(LC.res)
colnames(LC.res.f) <- paste0("LC.",colnames(LC.res.f))


#===========================================================================================
# merge expr 
#===========================================================================================
# the figure not suit set in Fig1 
#========================================
res.final <- cbind(LCBM.res.f,LC.res.f)
cor(res.final) -> cor.res
cor.res.f <- cor.res[grep("LCBM\\.",colnames(cor.res)),grep("LC\\.",colnames(cor.res))]









#===========================================================================================
# check CNV information with TCGA LUAD CNV level 
# just make a correlation plot 
#===========================================================================================

# first deal with TCGA LUAD info 
tcga.dat <- read.delim2("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/luad_tcga_pub_segments.seg",sep="\t",header=T)
tcga.dat$len <- tcga.dat$loc.end -tcga.dat$loc.start
tcga.dat$seg.mean <- as.numeric(as.vector(tcga.dat$seg.mean))
tcga.dat$score <- tcga.dat$len*tcga.dat$seg.mean

# need one more step to calculate values 
res.score <- aggregate(score~chrom,data=tcga.dat,FUN=sum)

res.len <- aggregate(len~chrom,data=tcga.dat,FUN=sum)


res.f <- merge(res.score,res.len,by="chrom")
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_CNV_seg.RDS")
# next deal with infercnv information 
res.ff <- res.f[-23,]






##############################################################################################
# 2021-2-24
# check CNV information with TCGA LUSC CNV level
#=============================================================================================

tcga.dat <- read.delim2("/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/lusc_tcga_segments.seg",sep="\t",header=T)
tcga.dat$len <- tcga.dat$loc.end -tcga.dat$loc.start
tcga.dat$seg.mean <- as.numeric(as.vector(tcga.dat$seg.mean))
tcga.dat$score <- tcga.dat$len*tcga.dat$seg.mean

# need one more step to calculate values 
res.score <- aggregate(score~chrom,data=tcga.dat,FUN=sum)

res.len <- aggregate(len~chrom,data=tcga.dat,FUN=sum)


res.f <- merge(res.score,res.len,by="chrom")
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/LUSC_CNV_seg.RDS")
# next deal with infercnv information 
res.ff <- res.f[-23,]



# Self scRNA data use scRNA to calculate ,code was in complexheatmap 
# you can find it by grep "corrplot" in complexheatmap.r
#==============================================================================================
# tmp.res 
tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/LUSC_CNV_seg.RDS")
tcga.res <- tcga.res[-23,]
tcga.res$value <- tcga.res$score/tcga.res$len
#===================================
# maybe the negative correlation sample has little tumor cells and little genes expressed
final.cor <- t(apply(tmp.res,2,function(x){
	tmpp <- cor.test(x-1,tcga.res$value)
	c(tmpp$p.value,tmpp$estimate)
}))

# pvalue is bigger than 0.05 ,and do not fit to plot ,just get a txt
final.cor.f <- final.cor[c(1:3),]
colnames(final.cor.f) <- c("pvalue","cor")
write.table(final.cor.f,file="/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_LUSC_BM_CNV_cor.txt",row.names=T,col.names=T,quote=F)





#============================================================================================
# 2020-12-22
# Use TCGA data to verify different stage have different CNV score 
#============================================================================================
options(stringsAsFactors=F)
tcga.dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/seg_based_scores.tsv",sep="\t",header=T)

# TCGA stage ann
tmp <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
ann <- tmp[,c(1,102:105)]
tmp.res <- merge(ann,tcga.dat,by.x="sampleID",by.y="Sample")
tmp.res$stage_type <- "unknow"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage I","Stage IA","Stage IB","Stage IIA","Stage IIB"))] <- "early"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))] <- "late"


# another try 
tmp.res$stage_type <- "unknow"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage I","Stage IA","Stage IB"))] <- "stageI"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage IIA","Stage IIB"))] <- "stageII"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage IIIA","Stage IIIB"))] <- "stageIII"
tmp.res$stage_type[which(tmp.res$pathologic_stage%in%c("Stage IV"))] <- "stageIV"
#result 
# Copy number burden scores frac_altered and n_segs (“fraction altered,” and “number of segments,” respectively) represent 
# the fraction of bases deviating from baseline ploidy (defined as above 0.1 or below −0.1 in log2 relative copy number (CN) space), 
# and the total number of segments in each sample’s copy number profile, respectively.
# So, lly think this score more high is more instability

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_Stage_CNV.pdf",useDingbats=F)
boxplot(frac_altered~stage_type,data=tmp.res,FUN=median)
legend("topright",
	legend=paste0("P= ",wilcox.test(tmp.res$frac_altered[which(tmp.res$stage_type=="late")],tmp.res$frac_altered[which(tmp.res$stage_type=="early")])$p.value))
dev.off()


#=========================================================================================
# another try 
# use this score to calculate correlation with BMS 
# 2020-12-22
#=========================================================
tcga.dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/seg_based_scores.tsv",sep="\t",header=T)
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
luad_mod <- data.frame(luad_mod)
luad_mod$ID <- gsub("\\.","-",rownames(luad_mod))














#===========================================================================================
# use GSE131907 to check 
#
#===========================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")
sub.dat <- subset(dat,cells=which(dat$Sample=="LUNG_T28"))
dat.f <- subset(sub.dat,cells=which(sub.dat$Cell_type.refined%in%c("Epithelial cells","T/NK cells")))

dir.create('/public/workspace/lily/Lung2Brain/GSE131907/Lung_T28')
#=========================
#===== just use llymarker
setwd("/public/workspace/lily/Lung2Brain/GSE131907/Lung_T28")


write.table(dat.f$Cell_type.refined,paste('Lung_T28',"cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
count<-as.matrix(dat.f@assays$RNA@counts)
write.table(count,paste("Lung_T28","count_exp.txt",sep='_'),sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste('Lung_T28',"count_exp.txt",sep='_'),
         annotations_file=paste('Lung_T28',"cell_info.txt",sep='_'),
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names=c("T/NK cells"))
         

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="./", 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             no_prelim_plot = TRUE,
                             num_threads=6, #big
                             no_plot=F ,
                             output_format = "pdf" # maybe can more quick 
                             # used for final scaling to fit range (0,2) centered at 1.
                             )

dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/Lung_T28/run.final.infercnv_obj")
Lung_T28 <- dat 
expr.Lung_T28 <- as.data.frame(Lung_T28@expr.data)
expr.Lung_T28$gene_order <- Lung_T28@gene_order$chr
res.Lung_T28 <- aggregate(.~gene_order,data=expr.Lung_T28,FUN=mean) #get every chr for each cell 
rownames(res.Lung_T28) <- res.Lung_T28$gene_order
res.Lung_T28$gene_order <- NULL
res.f.Lung_T28 <- apply(res.Lung_T28,1,function(x){mean(x)})

#==============================
tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_CNV_seg.RDS")










#================================================================================================
# tumor diversity 
# run each sample bu ROGUE 
# 2020-12-15
#================================================================================================
library(rlist)
library(Seurat)
library(pheatmap)
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))

#==================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
dat <- subset(dat,cells=which(dat$type=="maliganant"))
samplename <- unique(dat$orig.ident)
res.list <- list()
for(i in 1:length(samplename))
{
	tmp <- subset(dat,cells=which(dat$orig.ident==samplename[i]))
	expr <- as.matrix(tmp[['RNA']]@data)
#=============================================================================
# Run ROGUE 
#=============================================================================
	expr <- matr.filter(expr, min.cells = 10, min.genes = 10) #filter gene and cells 
	# random sampling for 10 times 
	tmp_res <- c()
	set.seed(12345)
	for(j in 1:10){
		subset.expr <- expr[,sample(1:ncol(expr),ncol(expr)/2)]
		ent.res <- SE_fun(subset.expr)
		rogue.value <- CalculateRogue(ent.res, platform = "UMI")
		tmp_res <- c(tmp_res,rogue.value)

	}
	res.list <- list.append(res.list,tmp_res)
	
}
names(res.list) <- samplename
saveRDS(res.list,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/Tumor_diversity.RDS")
lapply(res.list,function(x){1-x}) -> aa
# use each sample random 10 times result to boxplot 
c(aa[[1]],aa[[2]],aa[[7]]) -> LCBM
c(aa[[3]],aa[[4]],aa[[5]],aa[[6]]) -> LC

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Tumor_diversity.pdf")
boxplot(LCBM,LC,names=c("LCBM","LC"),outline=F)
legend("topright",legend=paste0("P=",wilcox.test(LCBM,LC)$p.value))
dev.off()








# 2020-12-17
#  line plot of CNV result 
#========================================================
# see in complexpheatmap.r





#========================================================
# LCBM 
# 2020-12-25
#========================================================
library(Seurat)
tmp <- readRDS("/public/workspace/lily/Lung2Brain/inte7/sigtumor_LC.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tmp[['RNA']]@data),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=0)









################################################################################################################################
# 2021-4-7
# use mouse lung development scRNA data to verfiy
#===============================================================================================================================
library(data.table)

E18 <- fread("~/metastasis/data/verify/GSE160876/GSM4885427_E18_JS-8_processed_counts.csv")
median(apply(E18[,-1],2,function(x){length(which(x[]>0))}))
P0_1 <- fread("~/metastasis/data/verify/GSE160876/GSM4885428_P0_JS-1_processed_counts.csv.gz")
P0_2 <- fread("~/metastasis/data/verify/GSE160876/GSM4885429_P0_JS-6_processed_counts.csv.gz")

P7_1 <- fread("~/metastasis/data/verify/GSE160876/GSM4885430_P7_JS-13_processed_counts.csv.gz")
P7_2 <- fread("~/metastasis/data/verify/GSE160876/GSM4885431_P7_JS-3_processed_counts.csv.gz")

P14_1 <- fread("~/metastasis/data/verify/GSE160876/GSM4885432_P14_JS-10_processed_counts.csv.gz")
P14_2 <- fread("~/metastasis/data/verify/GSE160876/GSM4885433_P14_JS-4_processed_counts.csv.gz")

P64_1 <- fread("~/metastasis/data/verify/GSE160876/GSM4885434_P64_JS-5_processed_counts.csv.gz")
P64_2 <- fread("~/metastasis/data/verify/GSE160876/GSM4885435_P64_JS-11_processed_counts.csv.gz")


median(apply(E18[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P0_1[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P0_2[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P7_1[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P7_2[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P14_1[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P14_2[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P64_1[,-1],2,function(x){length(which(x[]>0))}))
median(apply(P64_2[,-1],2,function(x){length(which(x[]>0))}))

# [1] 2636
# > median(apply(P0_1[,-1],2,function(x){length(which(x[]>0))}))
# [1] 1614
# > median(apply(P0_2[,-1],2,function(x){length(which(x[]>0))}))
# [1] 2159
# > median(apply(P7_1[,-1],2,function(x){length(which(x[]>0))}))
# [1] 2097
# > median(apply(P7_2[,-1],2,function(x){length(which(x[]>0))}))
# [1] 1571
# > median(apply(P14_1[,-1],2,function(x){length(which(x[]>0))}))
# [1] 1676
# > median(apply(P14_2[,-1],2,function(x){length(which(x[]>0))}))
# [1] 215
# > median(apply(P64_1[,-1],2,function(x){length(which(x[]>0))}))
# [1] 1113
# > median(apply(P64_2[,-1],2,function(x){length(which(x[]>0))}))
# [1] 1548







############################################################################################################################################
# 2021-4-7
# validate BRCA CNV info 
#===========================================================================================================================================
# tmp.res 
# prepare BRCA Mutation RDS
tcga.dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_BRCA/brca_tcga_pub2015_segments.seg",header=T)
tcga.dat$len <- tcga.dat$loc.end -tcga.dat$loc.start
tcga.dat$seg.mean <- as.numeric(as.vector(tcga.dat$seg.mean))
tcga.dat$score <- tcga.dat$len*tcga.dat$seg.mean

# need one more step to calculate values 
res.score <- aggregate(score~chrom,data=tcga.dat,FUN=sum)
res.len <- aggregate(len~chrom,data=tcga.dat,FUN=sum)
res.f <- merge(res.score,res.len,by="chrom")
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/TCGA_BRCA/BRCA_CNV_seg.RDS")


# 1. BRCA result 
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/inte7_chr_CNV.RDS")
tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_BRCA/BRCA_CNV_seg.RDS")
tcga.res <- tcga.res[-23,]
tcga.res$value <- tcga.res$score/tcga.res$len
#===================================
# maybe the negative correlation sample has little tumor cells and little genes expressed
final.cor <- t(apply(tmp.res,2,function(x){
	tmpp <- cor.test(x-1,tcga.res$value)
	c(tmpp$p.value,tmpp$estimate)
}))



#######################################################################################################################################
# 2. ME result  
# prepare Melanoma Mutation RDS
tcga.dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_SKCM/skcm_tcga_segments.seg",header=T)
tcga.dat$len <- tcga.dat$loc.end -tcga.dat$loc.start
tcga.dat$seg.mean <- as.numeric(as.vector(tcga.dat$seg.mean))
tcga.dat$score <- tcga.dat$len*tcga.dat$seg.mean

# need one more step to calculate values 
res.score <- aggregate(score~chrom,data=tcga.dat,FUN=sum)
res.len <- aggregate(len~chrom,data=tcga.dat,FUN=sum)
res.f <- merge(res.score,res.len,by="chrom")
saveRDS(res.f,file="/public/workspace/lily/metastasis/data/verify/TCGA_SKCM/SKCM_CNV_seg.RDS")


# 2. SKCM result 
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/inte7_chr_CNV.RDS")
tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_SKCM/SKCM_CNV_seg.RDS")
tcga.res <- tcga.res[-23,]
tcga.res$value <- tcga.res$score/tcga.res$len
#===================================
# maybe the negative correlation sample has little tumor cells and little genes expressed
final.cor <- t(apply(tmp.res,2,function(x){
	tmpp <- cor.test(x-1,tcga.res$value)
	c(tmpp$p.value,tmpp$estimate)
}))


#########################################################################################################################################
# 2021-4-7
# check TCGA LUAD CNV feature with TCGA BRCA CNV feature 
#========================================================================================================================================
tmp.res <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/inte7_chr_CNV.RDS")
tmp.dat <- tmp.res[,1:3]
tmp.dat <- tmp.dat-1

luad.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_CNV_seg.RDS")
luad.res <- luad.res[-23,]
luad.res$value <- luad.res$score/luad.res$len

skcm.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_SKCM/SKCM_CNV_seg.RDS")
skcm.res <- skcm.res[-23,]
skcm.res$value <- skcm.res$score/skcm.res$len


brca.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_BRCA/BRCA_CNV_seg.RDS")
brca.res <- brca.res[-23,]
brca.res$value <- brca.res$score/brca.res$len


# tmp.res 
lusc.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUSC/LUSC_CNV_seg.RDS")
lusc.res <- lusc.res[-23,]
lusc.res$value <- lusc.res$score/lusc.res$len



#===================================================================================================================================
# combine into a table

rs <- cbind(luad.res$value,skcm.res$value,lusc.res$value,brca.res$value)
colnames(rs) <- c("LUAD","SKCM","LUSC","BRCA")
rownames(rs) <- rownames(tmp.dat)



#===================================================================================================================================
# maybe use circlize to show correlation 
# 2021-4-16
####################################################################################################################################
a05 <- data.frame(t(apply(rs,2,function(x){
	tmp <- cor.test(x,tmp.dat$res.f.a05,method="spearman")
	c(tmp$estimate,tmp$p.value)
})))
a05$sample <- "A05"

a12 <- data.frame(t(apply(rs,2,function(x){
	tmp <- cor.test(x,tmp.dat$res.f.a12,method="spearman")
	c(tmp$estimate,tmp$p.value)
})))
a12$sample <- "A12"

tb <- data.frame(t(apply(rs,2,function(x){
	tmp <- cor.test(x,tmp.dat$res.f.tb,method="spearman")
	c(tmp$estimate,tmp$p.value)
})))
tb$sample <- "TB"

plot.rs <- rbind(a05,a12,tb)
rownames(plot.rs) <- NULL
colnames(plot.rs) <- c("Cor","Pvalue","Sample")
plot.rs$type <- rep(c("LUAD","SKCM","LUSC"),3)
plot.rs$logP <- -log2(plot.rs$Pvalue)
plot.rs$Sample <- factor(plot.rs$Sample,level=c("A12","A05","TB"))

# plot result by barplot 
library(ggplot2)

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Correlation_TCGA_CNV.pdf",useDingbats=F)
ggplot(plot.rs,aes(x=type,y=Cor,fill=logP,group=Sample)) + geom_bar(stat="identity",position="dodge") + theme_classic() + 
	geom_hline(aes(yintercept=0.5),colour="#990000") +
	scale_fill_gradientn(colours=c("#0091cd","white","#ecb731","#ed1b2e"))

dev.off()


































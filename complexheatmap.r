
#==================================================================
# change infercnv plot 
# use complexheatmap
#==================================================================
bytlib load R-3.6.0
bytlib load JAGS-4.3.0
R
library(Seurat)
library(infercnv)

#==================================================================
# try to use intersection gene to plot infercnv plot 
#==================================================================
dat.a05 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190305/run.final.infercnv_obj")
dat.a12 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190312/run.final.infercnv_obj")

# co_gene <- Reduce(intersect(rownames(dat.a05@expr.data),rownames(dat.a12@expr.data)))
co_gene <- intersect(rownames(dat.a05@expr.data),rownames(dat.a12@expr.data))
dat.merge <- cbind(dat.a05@expr.data[co_gene,],dat.a12@expr.data[co_gene,])
dat.res <- t(dat.merge)
# you have to check cell number 

library(ComplexHeatmap)
library(circlize)
# make a annotation 
# but maybe the gene position is not arrange 
gene_order <- dat.a05@gene_order[co_gene,]
col_chr <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B')
names(col_chr) <- paste0("chr",1:22)
bot_ann <- HeatmapAnnotation(df=gene_order[,1,drop=F],col=list(chr=col_chr))

col_sample <- c("blue","green")
names(col_sample) <- c("A20190305","A20190312")
tmp_cell <- data.frame(row.names=rownames(dat.res),sample=rep("unknow",nrow(dat.res)),stringsAsFactors=F)
# add some sample information
tmp_cell[which(rownames(tmp_cell)%in%colnames(dat.a05@expr.data)),"sample"] <- "A20190305"
tmp_cell[which(rownames(tmp_cell)%in%colnames(dat.a12@expr.data)),"sample"] <- "A20190312"


left_ann <- rowAnnotation(df=tmp_cell,col=list(samples=col_sample))



Heatmap(dat.res,col=colorRamp2(c(0,0.5,1,1.5,2),c("#00008B","#24249B","white","#9B2424","darkred")),
  name = "CNV", cluster_rows = FALSE, cluster_columns = FALSE,bottom_annotation = bot_ann,left_annotation=left_ann,
        show_row_names = FALSE, show_column_names = FALSE,border=TRUE,column_split=gene_order$chr,column_title=NULL,column_gap=unit(0,"mm")
)



#=========================================================================================
# infercnv score boxplot 
#=========================================================================================
library(infercnv)
#====================
dat.a05 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190305/run.final.infercnv_obj")
dat.a12 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190312/run.final.infercnv_obj")
dat.tb <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/T-Bsc1/run.final.infercnv_obj")
#
BT1296 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/BT1296/run.final.infercnv_obj")
BT1297 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/BT1297/run.final.infercnv_obj")
BT1431m <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/scrBT1431m/run.final.infercnv_obj")
BT1432m <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/scrBT1432m/run.final.infercnv_obj")





#===========
expr.a05 <- as.data.frame(dat.a05@expr.data[,-dat.a05@reference_grouped_cell_indices$T_cell])
expr.a05$gene_order <- dat.a05@gene_order$chr
res.a05 <- aggregate(.~gene_order,data=expr.a05,FUN=mean) #get every chr for each cell 
rownames(res.a05) <- res.a05$gene_order
res.a05$gene_order <- NULL
res.f.a05 <- apply(res.a05,1,function(x){mean(x)})

#===========
expr.a12 <- as.data.frame(dat.a12@expr.data[,-dat.a12@reference_grouped_cell_indices$T_cell])
expr.a12$gene_order <- dat.a12@gene_order$chr
res.a12 <- aggregate(.~gene_order,data=expr.a12,FUN=mean) #get every chr for each cell 
rownames(res.a12) <- res.a12$gene_order
res.a12$gene_order <- NULL
res.f.a12 <- apply(res.a12,1,function(x){mean(x)})

#===========
expr.tb <- as.data.frame(dat.tb@expr.data[,-dat.tb@reference_grouped_cell_indices$T_cell])
expr.tb$gene_order <- dat.tb@gene_order$chr
res.tb <- aggregate(.~gene_order,data=expr.tb,FUN=mean) #get every chr for each cell 
rownames(res.tb) <- res.tb$gene_order
res.tb$gene_order <- NULL
res.f.tb <- apply(res.tb,1,function(x){mean(x)})

#===========
expr.BT1296 <- as.data.frame(BT1296@expr.data[,-BT1296@reference_grouped_cell_indices$T_cell])
expr.BT1296$gene_order <- BT1296@gene_order$chr
res.BT1296 <- aggregate(.~gene_order,data=expr.BT1296,FUN=mean) #get every chr for each cell 
rownames(res.BT1296) <- res.BT1296$gene_order
res.BT1296$gene_order <- NULL
res.f.BT1296 <- apply(res.BT1296,1,function(x){mean(x)})

#===========
expr.BT1297 <- as.data.frame(BT1297@expr.data[,-BT1297@reference_grouped_cell_indices$T_cell])
expr.BT1297$gene_order <- BT1297@gene_order$chr
res.BT1297 <- aggregate(.~gene_order,data=expr.BT1297,FUN=mean) #get every chr for each cell 
rownames(res.BT1297) <- res.BT1297$gene_order
res.BT1297$gene_order <- NULL
res.f.BT1297 <- apply(res.BT1297,1,function(x){mean(x)})

#===========
expr.BT1431m <- as.data.frame(BT1431m@expr.data[,-BT1431m@reference_grouped_cell_indices$T_cell])
expr.BT1431m$gene_order <- BT1431m@gene_order$chr
res.BT1431m <- aggregate(.~gene_order,data=expr.BT1431m,FUN=mean) #get every chr for each cell 
rownames(res.BT1431m) <- res.BT1431m$gene_order
res.BT1431m$gene_order <- NULL
res.f.BT1431m <- apply(res.BT1431m,1,function(x){mean(x)})

#===========
expr.BT1432m <- as.data.frame(BT1432m@expr.data[,-BT1432m@reference_grouped_cell_indices$T_cell])
expr.BT1432m$gene_order <- BT1432m@gene_order$chr
res.BT1432m <- aggregate(.~gene_order,data=expr.BT1432m,FUN=mean) #get every chr for each cell 
rownames(res.BT1432m) <- res.BT1432m$gene_order
res.BT1432m$gene_order <- NULL
res.f.BT1432m <- apply(res.BT1432m,1,function(x){mean(x)})


#============================================================
# 
tmp.res <- cbind(data.frame(res.f.a05),data.frame(res.f.a12),data.frame(res.f.tb),data.frame(res.f.BT1296),data.frame(res.f.BT1297),
    data.frame(res.f.BT1431m),data.frame(res.f.BT1432m))

saveRDS(tmp.res,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/inte7_chr_CNV.RDS") # 2021-4-7 add sample cnv result 
# maybe can do some change to plot res 
library(reshape)
group.res <- apply(tmp.res[,1:7],1,function(x){c(median(x[1:3]),median(x[4:7]))})
rownames(group.res) <- c("LCBM","LC")
group.res.final <- melt(group.res,id="col.names")
colnames(group.res.final) <- c("group","chr","infercnv_score")
group.res.final$chr <- factor(group.res.final$chr,levels=unique(group.res.final$chr))
group.res.final$group <- factor(group.res.final$group,levels=unique(group.res.final$group))

pdf("/public/workspace/lily/Lung2Brain/inte7/infercnv/cnv_point.pdf",useDingbats=F)
ggplot(data = group.res.final,aes(x=chr,y=infercnv_score,group=group,color=group))+
	geom_line()+
	geom_point()+
	ylim(0.90,1.1)+
	theme_classic()+
	#theme(plot.title = element_text(hjust = 0.5))+
	theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1))
dev.off()

#============================================================
# compare with TCGA CNV info 

tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_CNV_seg.RDS")
tcga.res <- tcga.res[-23,]
tcga.res$value <- tcga.res$score/tcga.res$len
#===================================
# maybe the negative correlation sample has little tumor cells and little genes expressed
final.cor <- t(apply(tmp.res,2,function(x){
	tmpp <- cor.test(x-1,tcga.res$value)
	c(tmpp$p.value,tmpp$estimate)
}))

#===================================
# try to plot correlation figure 
#===================================
cor.tmp <- cbind(tcga.res[,"value"],tmp.res[,c("res.f.a05","res.f.a12","res.f.tb")]-1)
colnames(cor.tmp) <- c("TCGA_LUAD","A05","A12","TB_sc")
cor(cor.tmp) -> tmp

#colorRamp2(c(minvalue,(minvalue+1)/2,1,(maxvalue+1)/2,maxvalue),c("#00008B","#24249B","white","#9B2424","darkred"))
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TCGA_BM_cor.pdf")
corrplot(tmp, method = "circle", type = "upper", tl.pos = "d",cl.lim=c(0,1))
corrplot(tmp, add = TRUE, type = "lower", method = "number", diag = FALSE, tl.pos = "n", cl.pos = "n",cl.lim=c(0,1))
dev.off()

#=====================================================================================================================
# another way to show cnv inforamtion 
# use mad to get boxplot 
#=====================================================================================================================
# library(reshape)
# tmp.res$chr <- rownames(tmp.res)
# res.final <- melt(tmp.res,id.vars="chr")
# res.final$group <- "unknow"
# colnames(res.final) <- c("chr","sample","infercnv_score","group")
# res.final$group[which(res.final$sample%in%c("res.f.a05","res.f.a12","res.f.tb"))] <- "LCBM"
# res.final$group[which(res.final$sample%in%c("res.f.BT1296","res.f.BT1297","res.f.BT1431m","res.f.BT1432m"))] <- "LC"

# #======================================================
# #======================================================
# library(ggplot2)

# # some prepare for ggplot
# res.final$chr <- factor(res.final$chr,levels=unique(res.final$chr))
# res.final$group <- factor(res.final$group,levels=unique(res.final$group))

# pdf("/public/workspace/lily/Lung2Brain/inte7/infercnv/infercnv_score.pdf",useDingbats=F)
# ggplot(data = res.final,aes(x=chr,y=infercnv_score,fill=group))+ geom_boxplot()+
#     geom_boxplot(outlier.size=0)+
# 	ylim(0.9,1.2)+
# 	theme_classic()+
# 	theme(plot.title = element_text(hjust = 0.5))
# dev.off()
#==================================================================================================================
# use single cell data to calculate 
mad.res.a05 <- apply(res.a05,2,function(x){mad(x)})
mad.res.a12 <- apply(res.a12,2,function(x){mad(x)})
mad.res.tb <- apply(res.tb,2,function(x){mad(x)})
mad.res.BT1296 <- apply(res.BT1296,2,function(x){mad(x)})
mad.res.BT1297 <- apply(res.BT1297,2,function(x){mad(x)})
mad.res.BT1431m <- apply(res.BT1431m,2,function(x){mad(x)})
mad.res.BT1432m <- apply(res.BT1432m,2,function(x){mad(x)})

pdf("/public/workspace/lily/Lung2Brain/inte7/infercnv/MAD_cnv_plot.pdf",useDingbats=F)
boxplot(c(mad.res.a05,mad.res.a12,mad.res.tb),c(mad.res.BT1296,mad.res.BT1297,mad.res.BT1431m,mad.res.BT1432m),
    names=c("LCBM","LC"),ylabs="MAD",outline=F)
legend('topright',legend=paste0("P= ",round(wilcox.test(c(mad.res.a05,mad.res.a12,mad.res.tb),c(mad.res.BT1296,mad.res.BT1297,mad.res.BT1431m,mad.res.BT1432m))$p.value,digits=4)))
dev.off()















#========================================================================================================
# use complex heatmap to plot result 
# plot each chr and then cbind the result 
#========================================================================================================
bytlib load R-3.6.0
bytlib load JAGS-4.3.0
R
library(Seurat)
library(infercnv)
library(ComplexHeatmap)
library(circlize)
options(stringsAsFactors=F)
#====================
dat.a05 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190305/run.final.infercnv_obj")
dat.a12 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190312/run.final.infercnv_obj")
dat.tb <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/T-Bsc1/run.final.infercnv_obj")
#
BT1296 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/BT1296/run.final.infercnv_obj")
BT1297 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/BT1297/run.final.infercnv_obj")
BT1431m <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/scrBT1431m/run.final.infercnv_obj")
BT1432m <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/scrBT1432m/run.final.infercnv_obj")


#=====================
# chr 1 
# test 
#=====================
#------------------
expr.a05 <- as.data.frame(dat.a05@expr.data)
#expr.a05$gene_order <- dat.a05@gene_order$chr
#------------------
expr.a12 <- as.data.frame(dat.a12@expr.data)
#expr.a12$gene_order <- dat.a12@gene_order$chr
#------------------
expr.tb <- as.data.frame(dat.tb@expr.data)
#expr.tb$gene_order <- dat.tb@gene_order$chr
#------------------
expr.BT1296 <- as.data.frame(BT1296@expr.data)
#expr.BT1296$gene_order <- BT1296@gene_order$chr
#------------------
expr.BT1297 <- as.data.frame(BT1297@expr.data)
#expr.BT1297$gene_order <- BT1297@gene_order$chr
#------------------
expr.BT1431m <- as.data.frame(BT1431m@expr.data)
#expr.BT1431m$gene_order <- BT1431m@gene_order$chr
#------------------
expr.BT1432m <- as.data.frame(BT1432m@expr.data)
#expr.BT1432m$gene_order <- BT1432m@gene_order$chr

expr.list <- list("expr.a05" = as.data.frame(dat.a05@expr.data[,-dat.a05@reference_grouped_cell_indices$T_cell]),
		"expr.a12" = as.data.frame(dat.a12@expr.data[,-dat.a12@reference_grouped_cell_indices$T_cell]),
		"expr.tb" = as.data.frame(dat.tb@expr.data[,-dat.tb@reference_grouped_cell_indices$T_cell]),
		"expr.BT1296" = as.data.frame(BT1296@expr.data[,-BT1296@reference_grouped_cell_indices$T_cell]),
		"expr.BT1297" = as.data.frame(BT1297@expr.data[,-BT1297@reference_grouped_cell_indices$T_cell]),
		"expr.BT1431m" = as.data.frame(BT1431m@expr.data[,-BT1431m@reference_grouped_cell_indices$T_cell]),
		"expr.BT1432m" = as.data.frame(BT1432m@expr.data[,-BT1432m@reference_grouped_cell_indices$T_cell])
	)

# get_rownames <- function(x,y){
# 	if( !is.null(rownames(x)) ){
# 		res <- union(rownames(x),rownames(y))
# 	}else{
# 		res <- union(x,rownames(y))
# 	}
	
# 	return(res)
# }

all_join <- function(x,y){
	res <- merge(x,y,by="row.names",all.x=T,all.y=T)
	res$Row.names  -> rownames(res)
	res$Row.names <- NULL
	return(res)

}

# all_genes <- purrr::reduce(expr.list,.f=get_rownames)

# for(i in 1:length(expr.list)){
# 	tmp <- expr.list[[i]]
# 	genes <- all_genes[-which(all_genes%in%rownames(tmp))]

# }

res.f <- purrr::reduce(expr.list,all_join)
# maybe not ok 
# res.f[which(is.na(res.f))] <- 1
#=========================================================
res.final <- apply(res.f,2,function(x){
	tmp <- x
	tmp[which(is.na(tmp))] <- 1
	c(tmp)
})

saveRDS(res.final,file="/public/workspace/lily/Lung2Brain/inte7/infercnv/all_cell_mat.RDS")
saveRDS(expr.list,file="/public/workspace/lily/Lung2Brain/inte7/infercnv/expr.list.RDS")
#==========================================================
# check infercnv with TCGA 
#==========================================================
chr.info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",sep="\t")
colnames(chr.info) <- c("gene","chr","gene_start","gene_end")

tmp.res <- sapply(expr.list,function(x){
	tmp <- x
	tmp$gene <- rownames(tmp)
	res.f <- merge(tmp,chr.info[,1:2],by="gene")
	# trans gene names to row names
	rownames(res.f) <- res.f$gene
	res.f$gene <- NULL
	res.ff <- aggregate(.~chr,data=res.f,FUN=mean) #get every chr for each cell 
	rownames(res.ff) <- res.ff$chr
	res.ff$chr <- NULL
	res.final <- apply(res.ff,1,function(x){mean(x)})
	c(res.final)
})

tmp.res <- data.frame(tmp.res)
tmp.res$chr <- as.numeric(gsub("^chr","",rownames(tmp.res)))
tmp.res <- tmp.res[order(tmp.res$chr),]


#==========
# load tcga cnv info
tcga.res <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_CNV_seg.RDS")
tcga.res <- tcga.res[-23,]

#===================================
# maybe the negative correlation sample has little tumor cells and little genes expressed
final.cor <- t(apply(tmp.res[,-8],2,function(x){
	tmpp <- cor.test(x-1,tcga.res$value)
	c(tmpp$p.value,tmpp$estimate)
}))







#=====================================================================================================================
# ready to plot result
#==============================================================
library(ComplexHeatmap)
library(circlize)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/all_cell_mat.RDS")
dat <- t(dat) 	# should be row is cells and col is gene , so should be 10321 cells 
# expr.list <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/expr.list.RDS")
# get gene and chromsome info 
chr.info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",sep="\t")
colnames(chr.info) <- c("gene","chr","gene_start","gene_end")

# # get cell for each samples 
# #==============================================================
# library(Seurat)
# #==============================================================
# library(ComplexHeatmap)
# #plot chr1   
# tmp.chr1 <- dat[which(rownames(dat)%in%chr.info$gene[which(chr.info$chr=="chr1")]),]
# tmp.chr1.f <- tmp.chr1-1

#===================================================================================================================
# make a annotation 

col_chr <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B')
names(col_chr) <- paste0("chr",1:22)

#=============================================
# add some row annotation 
#=============================================
col_sample <- c("#0081a7","#00afb9","#fdfcdc","#fed9b7","#f07167","#98c1d9","#d4e09b")
names(col_sample) <- c("A20190305","A20190312","T_Bsc","BT1296","BT1297","BT1431m","BT1432m")
tmp_cell <- data.frame(row.names=rownames(dat),sample=rep("unknow",nrow(dat)),stringsAsFactors=F)
# set cells names for each sample 
expr.list <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/expr.list.RDS")
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.a05)),"sample"] <- "A20190305"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.a12)),"sample"] <- "A20190312"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.tb)),"sample"] <- "T_Bsc"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.BT1296)),"sample"] <- "BT1296"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.BT1297)),"sample"] <- "BT1297"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.BT1431m)),"sample"] <- "BT1431m"
tmp_cell[which(rownames(tmp_cell)%in%colnames(expr.list$expr.BT1432m)),"sample"] <- "BT1432m"
left_ann <- rowAnnotation(df=tmp_cell,col=list(samples=col_sample))


#=============================================
# set some col information
#=============================================
tmp.chr <- chr.info[which(chr.info$gene%in%colnames(dat)),"chr",drop=F]
rownames(tmp.chr) <- NULL
bot_ann <- HeatmapAnnotation(df=tmp.chr[,1,drop=F],col=list(chr=col_chr)) # set each chromsome colors


#=================================================================================================
# plot result 
#=================================================================================================
pdf("/public/workspace/lily/Lung2Brain/inte7/infercnv/gene_union.pdf")
Heatmap(dat,col=colorRamp2(c(0,0.5,1,1.5,2),c("#00008B","#24249B","white","#9B2424","darkred")),
  name = "CNV", cluster_rows = FALSE, cluster_columns = FALSE,bottom_annotation = bot_ann,left_annotation=left_ann,
        show_row_names = FALSE, show_column_names = FALSE,border=TRUE,column_split=tmp.chr$chr,column_title=NULL,column_gap=unit(0,"mm")
)
dev.off()













































#================================================================================================
# try to another plot 
# plot each sample and each chr is same length
# plot each sample and then paste 
#================================================================================================
width<-c(13.97,13.688,11.148,10.865,10.16,9.596,9.031,8.325,7.903,7.62,7.62,7.62,6.491,6.067,5.786,5.08,4.516,4.374,3.387,3.528,2.681,2.963)

# you have to plot each result for each chr 
expr.list <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/expr.list.RDS")
chr.info <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",sep="\t")
colnames(chr.info) <- c("gene","chr","gene_start","gene_end")
# color for each chr 
# col_chr <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
#                   '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B')
# names(col_chr) <- paste0("chr",1:22)

# sample BT1296

for(i in 1:length(expr.list)){

	tmp.dat <- t(expr.list[[i]])
	tmp.chr <- chr.info[which(chr.info$gene%in%colnames(tmp.dat)),c("chr","gene"),drop=F]
	rownames(tmp.chr) <- NULL
	# set a more beatuiful color 
	ht_opt(heatmap_border = TRUE)
	ht_list <- list()
	for(j in 1:22){
		# subset data 
		gene.chr <- as.vector(tmp.chr$gene[which(tmp.chr$chr==paste0("chr",j))])
		hpdat <- tmp.dat[,gene.chr]-1
		anno_col<-data.frame(row.names=gene.chr,chr=rep(paste('chr',j,sep=''),length(gene.chr)),stringsAsFactors=F)
		cl<-c('white','black')
		se=j%%2+1
		col_chr<-cl[se]
		names(col_chr)<-paste('chr',j,sep='')
		# set ann col 
		tmp.bot_ann <- HeatmapAnnotation(df=anno_col,col=list(chr=col_chr),show_legend=FALSE,show_annotation_name=FALSE,border=TRUE)
		
		ht <- Heatmap(hpdat,col=colorRamp2(c(min(hpdat)-0.1,min(hpdat)/2,0,max(hpdat)/2,max(hpdat)+0.1),c("#00008B","#24249B","white","#9B2424","darkred")),use_raster=T,bottom_annotation=tmp.bot_ann,
				name = "CNV", cluster_rows = FALSE, cluster_columns = FALSE,width=unit(6,"mm"),
				show_row_names = FALSE, show_column_names = FALSE,column_title=j,column_title_side="bottom",border="#4d4f53")
		ht_list[[paste("ht",j,sep="")]]<-ht
	}
	ht_list.f <- Reduce("+", ht_list)
	pdf(paste0("/public/workspace/lily/Lung2Brain/inte7/infercnv/",names(expr.list)[i],"_infercnv.pdf"),width=15,height=12)
	draw(ht_list.f, ht_gap = unit(0,"mm"))
	dev.off()
}


ht_list <- list()  
ht_opt(heatmap_border = TRUE)
for (i in 1:22) 
{
 gene<-rownames(infercnv_obj@gene_order)[which(infercnv_obj@gene_order$chr==(paste('chr',i,sep='')))]
 hpdat<-tmp2[,gene]
 anno_col<-data.frame(row.names=gene,chr=rep(paste('chr',i,sep=''),length(gene)),stringsAsFactors=F)
 cl<-c('white','black')
 se=i%%2+1
 col_chr<-cl[se]
 names(col_chr)<-paste('chr',i,sep='')
 ha_col = HeatmapAnnotation(df=anno_col,col=list(chr=col_chr),show_legend=FALSE,show_annotation_name=FALSE,border=TRUE)

    ht <- Heatmap(hpdat,col=colorRamp2(c(-1,-0.5,0,0.5,1),c("#00008B","#24249B","white","#9B2424","darkred")),
            name = "CNV", cluster_rows = FALSE, cluster_columns = FALSE,width=unit(width[i],"mm"),
   bottom_annotation = ha_col,
            show_row_names = FALSE, show_column_names = FALSE,column_title=i,column_title_side="bottom",border=TRUE)
    ht_list[[paste("ht",i,sep="")]]<-ht
}
ht_list2 <- Reduce("+", ht_list)

pdf("infervnc_20201112.pdf",width=10)
draw(ht_list2, ht_gap = unit(0,"mm"))
dev.off()


#=======================================================================================================================
# WES result 
# plot =========================================
# A05
dat.05 <- read.table('/public/workspace/lily/Lung2Brain/WES/Sequenza/A20190305/A20190305_segments.txt',sep="\t",header=T,stringsAsFactors=F)
dat.12 <- read.table('/public/workspace/lily/Lung2Brain/WES/Sequenza/A20190312/A20190312_segments.txt',sep="\t",header=T,stringsAsFactors=F)
dat.tb <- read.table('/public/workspace/lily/Lung2Brain/WES/Sequenza/T_Bsc/T_Bsc_segments.txt',sep="\t",header=T,stringsAsFactors=F)

































#########################################################################################################################################################
# 2021-4-10
# check 8q genes 
#========================================================================================================================================================
library(infercnv)
# load GRch37 data

load("/public/workspace/lily/REF/CpRMAP.RData")

# sample A20190305
dat.a05 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190305/run.final.infercnv_obj")
expr.a05 <- as.data.frame(dat.a05@expr.data[,-dat.a05@reference_grouped_cell_indices$T_cell])
tmp.a05 <- apply(expr.a05,1,function(x){mean(x)})
tmp.res.a05 <- data.frame(gene=names(tmp.a05),value=unname(tmp.a05))
res.a05 <- merge(GRCh37.chrpqb[,c("GENE","CHR")],tmp.res.a05,by.x="GENE",by.y="gene")


# sample A20190312
dat.a12 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/A20190312/run.final.infercnv_obj")
expr.a12 <- as.data.frame(dat.a12@expr.data[,-dat.a12@reference_grouped_cell_indices$T_cell])
tmp.a12 <- apply(expr.a12,1,function(x){mean(x)})
tmp.res.a12 <- data.frame(gene=names(tmp.a12),value=unname(tmp.a12))
res.a12 <- merge(GRCh37.chrpqb[,c("GENE","CHR")],tmp.res.a12,by.x="GENE",by.y="gene")


# sample T_Bsc
dat.tb <- readRDS("/public/workspace/lily/Lung2Brain/inte7/infercnv/T-Bsc1/run.final.infercnv_obj")
expr.tb <- as.data.frame(dat.tb@expr.data[,-dat.tb@reference_grouped_cell_indices$T_cell])
tmp.tb <- apply(expr.tb,1,function(x){mean(x)})
tmp.res.tb <- data.frame(gene=names(tmp.tb),value=unname(tmp.tb))
res.tb <- merge(GRCh37.chrpqb[,c("GENE","CHR")],tmp.res.tb,by.x="GENE",by.y="gene")



###############################################################################################################################
tmp.1 <- res.a05[grep("8q",res.a05$CHR),]
tmp.2 <- res.a12[grep("8q",res.a12$CHR),]
tmp.3 <- res.tb[grep("8q",res.tb$CHR),]

tmp.1.f <- tmp.1[which(tmp.1$value>1),]
tmp.2.f <- tmp.2[which(tmp.2$value>1),]
tmp.3.f <-tmp.3[which(tmp.3$value>1),]


gene <- intersect(intersect(tmp.1.f$GENE,tmp.2.f$GENE),tmp.3.f$GENE)
# write.table(gene,file="./tmp.txt",quote=F)
# this signature is not ok 

tmp.res <- cbind(tmp.1.f[which(tmp.1.f$GENE%in%gene),],tmp.2.f[which(tmp.2.f$GENE%in%gene),],tmp.3.f[which(tmp.3.f$GENE%in%gene),])
tmp.res$value.mean <- apply(tmp.res[,c(3,6,9)],1,function(x){median(x)})
tmp.res <- tmp.res[order(tmp.res$value.mean,decreasing=T),]  # use value > 1.10 genes 
EIF3E
OXR1
SLC25A32
DCAF13
LRP12
ATP6V1C1
EMC2
NUDCD1
AZIN1
# this gene signature is OK , OS is significant and PFS is not signifcant but have tendency.
# and we also check the MYC gene expression ,found this gene is up reguation in three sample 

#=========================================================================================================================================
dat.05 <- read.table('/public/workspace/lily/Lung2Brain/WES/Varscan/A20190305/A20190305.seg',sep="\t",header=T,stringsAsFactors=F)
dat.12 <- read.table('/public/workspace/lily/Lung2Brain/WES/Varscan/A20190312/A20190312.seg',sep="\t",header=T,stringsAsFactors=F)
dat.tb <- read.table('/public/workspace/lily/Lung2Brain/WES/Varscan/T_Bsc/T_Bsc.seg',sep="\t",header=T,stringsAsFactors=F)

# summary(dat.05$log2_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -7.8330 -0.5510 -0.2750 -0.3224 -0.0430  2.7490
# -0.55 -0.1
# summary(dat.12$log2_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# -4.5890  0.1700  0.3910  0.3494  0.5820  3.8740
# 0.17 0.58
# summary(dat.tb$log2_ratio)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# -4.22000 -0.30500  0.06300  0.02168  0.38100  4.35600
#-0.3 0.38













































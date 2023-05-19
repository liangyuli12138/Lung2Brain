#!/usr/bin/Rscript
# 2021-6-28
# this program is used to add some plot for Figure1 
# the main body for figure1 is in Inte7_Figure1.R 


# 0. LUAD gene signature https://www.nature.com/articles/1209442
# come from nature medicine
library(Seurat)
dat <-readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
# tumor <- subset(dat,cells=which(dat$type=="maliganant"))

# calculate 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("LUAD_tumor","Lung_gene"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version5/Fig1/LUAD_sig_mod.RDS")


# plot result 
mod$group <- dat$maliganant
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig1/Res/LUAD_tumor_sig.pdf",useDingbats=F)
boxplot(LUAD_tumor_norm~group,data=mod,FUN=median,outline=F,ylims=c(0,1))
dev.off()








# replot one result 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
data <- dat[['RNA']]@data

tumor_gene <- apply(data[,which(dat$maliganant=="tumor")],2,function(x){length(which(x[]>0))})
ntumor_gene <- apply(data[,which(dat$maliganant=="non-tumor")],2,function(x){length(which(x[]>0))})

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig1/Res/tumor_ntumor_gene.pdf")
boxplot(tumor_gene,ntumor_gene,names=c('malignant','non-malignant'),main="express gene numbers",col=c("#E41A1C","#A6CEE3"),outline=F,ylim=c(0,8000))
legend('topright',legend=paste('p =',round(wilcox.test(tumor_gene,ntumor_gene)$p.value,5)),bty='n')
dev.off()












# 0. change a idea ,use all cells to calculate inferCNV 
#================================================================================================================================================
library(infercnv)
#================================================================================================================================================

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
dat$maliganant[which(dat$maliganant=="non-tumor")] <- "non_tumor"
sample <- unique(dat$orig.ident)

for(i in 1:length(sample)){
	tmp <- subset(dat,cells=which(dat$orig.ident==sample[i]))
	dir.create(paste0("/public/workspace/lily/Lung2Brain/Version5/infercnv/",sample[i]))
	setwd(paste0("/public/workspace/lily/Lung2Brain/Version5/infercnv/",sample[i]))

	write.table(tmp$maliganant,paste(sample[i],"cell_info.txt",sep='_'),sep="\t",col.names=F,quote=F)
	count<-as.matrix(tmp@assays$RNA@counts)
	write.table(count,paste(sample[i],"count_exp.txt",sep='_'),sep="\t",quote=F)

	infercnv_obj = CreateInfercnvObject(raw_counts_matrix=paste(sample[i],"count_exp.txt",sep='_'),
         annotations_file=paste(sample[i],"cell_info.txt",sep='_'),
         delim="\t",
         gene_order_file="/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt",
         ref_group_names=c("non_tumor"))
         

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


}


# check result
# get sample expr data
dat.list <- list()
file <- dir("/public/workspace/lily/Lung2Brain/Version5/infercnv/")
for(i in 1:length(file)){
    tmp <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Version5/infercnv/",file[i],"/run.final.infercnv_obj"))
    tmp.exp <- tmp@expr.data
    tmp.a <- apply(tmp.exp,2,function(x){sqrt(sum((x-1)^2))})
    tmp.rs <- data.frame(row.names=names(tmp.a),CNV.score=unname(tmp.a),sample=rep(file[i],length(tmp.a)))
    dat.list[[i]] <- tmp.rs
    names(dat.list)[i] <- file[i]
}

# then combine all cells 
rs <- dat.list[[1]]
rs <- rbind(rs,dat.list[[2]])
rs <- rbind(rs,dat.list[[3]])
rs <- rbind(rs,dat.list[[4]])
rs <- rbind(rs,dat.list[[5]])
rs <- rbind(rs,dat.list[[6]])
rs <- rbind(rs,dat.list[[7]])
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
rs <- rs[colnames(dat),]
rs$type <- dat$type
saveRDS(rs,file="/public/workspace/lily/Lung2Brain/Version5/Fig1/inte7_CNVscore.RDS")

# get Seurat object andplot result 
rs <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig1/inte7_CNVscore.RDS")
aggregate(CNV.score~type,data=rs,FUN=median)

# Plot
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(gg.gap)

rs$type <- factor(rs$type,levels=c("maliganant","Oligodendrocyte","Fibroblast","Endothelial","Myeloid","T_cell","B_cell"))
p1 <- ggplot(data=rs,aes(x=type, y=CNV.score, fill=type)) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.2, color="black", alpha=0.2,outlier.shape=NA) +
    scale_fill_viridis(discrete = TRUE) +
    theme( legend.position="none") + ylim(0,10)
    

p2 <- gg.gap(plot = p1,
    segments = c(10, 20),
    tick_width = 10, # 坐标间隔
    rel_heights = c(0.5,0, 0.25),# 设置分隔为的三个部分的宽度
    ylim = c(0,30)
)

pdf("/public/workspace/lily/Lung2Brain/Version5/Fig1/Res/CNV.score.celltype.pdf",useDingbats=F)
print(p1)

ggplot(data=rs,aes(x=type, y=CNV.score, fill=type)) +
    geom_violin(width=1.4) +
    geom_boxplot(width=0.2, color="black", alpha=0.2,outlier.shape=NA) +
    scale_fill_viridis(discrete = TRUE) +
    theme( legend.position="none") 
    
dev.off()














# 1. calculate Tumor between LCBM and LC for cytotrace 
# 2021 -6- 29
library(Seurat)
load("~/Lung2Brain/inte7/res_Data/cytotrace_result.RData")
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")

all(names(results$CytoTRACE)==colnames(tumor))
pdf("/public/workspace/lily/Lung2Brain/Version5/Fig1/Res/Cytotrace_Tumor.pdf",useDingbats=F)
boxplot(CytoTrace~type_group,data=tumor@meta.data,FUN=median,outline=F,main="Stemness")
dev.off()





# 2. Caculate entropy 
# not very good ,myedloid show highest ocuupy
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")

# calculate occpuicy
clusters <- as.numeric(as.vector(unique(dat$seurat_clusters)))
dat$orig.ident <- factor(dat$orig.ident,levels=unique(dat$orig.ident))
tmp <- sapply(clusters,function(i){
    
    idx <- unname(which.max(table(dat$orig.ident[which(dat$seurat_clusters==i)])/table(dat$orig.ident)))
    rs <- unname(table(dat$orig.ident[which(dat$seurat_clusters==i)])[idx])/length(which(dat$seurat_clusters==i))
    return(rs)
})

tmp.res <- data.frame(occupy= tmp,cluster= clusters)

tmp.res$type <- "unknow"
tmp.res$type[which(tmp.res$cluster%in%c(1,2,3,18,39,36,9,30,19,34,40,43,15,28,35,38))] <- "T_cell"
tmp.res$type[which(tmp.res$cluster%in%c(4,6,12,13,16,7,26,25,14))] <- "Myeloid"
tmp.res$type[which(tmp.res$cluster%in%c(11,23))] <- "B_cell"
tmp.res$type[which(tmp.res$cluster%in%c(29,27))] <- "Oligodendrocyte"
tmp.res$type[which(tmp.res$cluster%in%c(22,33))] <- "Fibroblast"
tmp.res$type[which(tmp.res$cluster%in%c(32))] <- "Endothelial"
tmp.res$type[which(tmp.res$cluster%in%c(0,5,10,20,21,37,24,17,8,31))] <- "malignant"



plot(density(tmp.res[which(tmp.res$type=="malignant"),"occupy"]),col="red")
lines(density(tmp.res[which(tmp.res$type=="T_cell"),"occupy"]),col="green")
lines(density(tmp.res[which(tmp.res$type=="Myeloid"),"occupy"]),col="orange")










#============================================================================================================================================
# TCGA LUAD chr8 verify 
# 2021-12-17

tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Data/luad_tcga_pub_segments.seg",header=T,sep="\t")

dat <- tmp.dat[which(tmp.dat$chrom==8&tmp.dat$loc.start>46000123),]
dat$length <- dat$loc.end - dat$loc.start
dat$value <- dat$seg.mean * dat$length
tmp.res <- aggregate(.~ID,data=dat[,c(1,7,8)],FUN=sum)
tmp.res$res <- tmp.res$value / tmp.res$length
tmp.res$ID <- gsub("-",".",tmp.res$ID)
tmp.res$status <- ifelse(tmp.res$res>0,"gain","loss")
# load BMS update score 
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/TCGA_LUAD_mod.RDS")
mod$ID <- rownames(mod)
res <- merge(tmp.res,mod[,c(2,4)],by="ID")
#



































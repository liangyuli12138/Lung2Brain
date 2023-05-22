#!/usr/bin/Rscript
################################### Metabolic#########################
library(Seurat)
library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(ggrepel)
options(stringsAsFactors=FALSE)
library(gtools)
library(scran)

source('/public/workspace/lily/software/SingleCellMetabolic/utils.R')
source('/public/workspace/lily/software/SingleCellMetabolic/runGSEA_preRank.R')


# read data
# all cell is tumor
tmp_data=readRDS("/public/workspace/lily/Lung2Brain/Version5/Data/Tumor_pur/pur_LCBM_clust.RDS")
all_data=as.matrix(tmp_data[['RNA']]@data) # make sure your data is correct size which can use as.matrix
cell_type=as.vector(tmp_data$seurat_clusters)  # cell type which means cluster or cell type 
tumor=unname(tmp_data$orig.ident) # this do not need to change ,this means all cells need to analysis



#=======================================================================================================================
col_data <- data.frame(tumor=tumor,cellType=as.character(cell_type),row.names=names(cell_type))
pathways <- gmtPathways("/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE
#set SingleCellExperiment object which include expr matrix ,cluster info ,and gene info
sce <- SingleCellExperiment(
  assays = all_data,
  colData = col_data,
  rowData = row_data
)

selected_tumor_sce <- sce #the example code use "selected_tumor_sce" as the name
selected_tumor_metabolic_sce <- sce[rowData(sce)$metabolic,] # dims :1506 37637; 1506 metabolic genes

######pheatmap show metabolic sce object
# ###### can not alloate size ?
# dat <- assay(selected_tumor_metabolic_sce)
# dist_dat  <- cor(dat,method="spearman")
# hc <- hclust(as.dist(1-dist_dat),method="ward.D2")
# mycolor = colorRampPalette(c("white","red"))(10)
# mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)
# col_row <- as.data.frame(colData(selected_tumor_metabolic_sce)[,"tumor",drop=F])
# col_row$tumor <- factor(col_row$tumor)
# pdf("/public/workspace/lily/res_cor_inte/Metabolic/malignant_metabolic_correlationMatrix.pdf",width=5,height=4,onefile=T)
# pheatmap(dist_dat,show_rownames = F,show_colnames = F,cluster_rows = hc,cluster_cols =hc,annotation_row = col_row,color=mycolor, annotation_legend = T)
# dev.off()

# metabolic landscape
# set.seed(12345)
# tsne_out <- Rtsne(t(assay(selected_tumor_metabolic_sce)),initial_dims=20,theta=0.0,perplexity = 30)
# tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_tumor_metabolic_sce)$tumor)
# tmp$cluster=factor(gsub('[0-9]*','',tmp$group)) # delete the number of sample names
# g <- ggplot(tmp) + geom_point(aes(x, y, colour = cluster), size = 1) +
#   labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
# ggsave("/public/workspace/lily/res_cor_inte/Metabolic/tumor_metabolic_tsne_PM_PwoM_LM.png",g,width=4,height=3)




###################################### scRNA_pathway_activity ##########################################
pathway_file <- "/public/workspace/lily/software/SingleCellMetabolic/Data/KEGG_metabolism.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
all_cell_types <- as.vector(selected_tumor_metabolic_sce$cellType)
cell_types <- unique(all_cell_types)

gene_pathway_number <- num_of_pathways(pathway_file,rownames(selected_tumor_metabolic_sce)[rowData(selected_tumor_metabolic_sce)$metabolic])
set.seed(123)
normalization_method <- "Deconvolution"

##Calculate the pathway activities
#mean ratio of genes in each pathway for each cell type
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = list(pathway_names,cell_types))
###calculate the pvalues using shuffle method
pvalues_mat <- matrix(NA,nrow=length(pathway_names),ncol=length(cell_types),dimnames = (list(pathway_names, cell_types)))

norm_tpm <- all_data


for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(norm_tpm))
  if(length(genes_comm) < 5) next

  pathway_metabolic_tpm <- norm_tpm[genes_comm, ]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[rowSums(pathway_metabolic_tpm)>0,]
  
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))

  #remove genes which are zeros in any celltype to avoid extreme ratio value
  keep <- colnames(mean_exp_eachCellType)[colAlls(mean_exp_eachCellType>0.001)]

  if(length(keep)<3) next
  
  #using the loweset value to replace zeros for avoiding extreme ratio value
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_metabolic_tpm <- t( apply(pathway_metabolic_tpm,1,function(x) {x[x<=0] <- min(x[x>0]);x} ))

  
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  #
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  #exclude the extreme ratios
  col_quantile <- apply(ratio_exp_eachCellType,2,function(x) quantile(x,na.rm=T))
  col_q1 <- col_quantile["25%",]
  col_q3 <- col_quantile["75%",]
  col_upper <- col_q3 * 3
  col_lower <- col_q1 / 3
  outliers <- apply(ratio_exp_eachCellType,1,function(x) {any( (x>col_upper)|(x<col_lower) )} )
  
  if(sum(!outliers) < 3) next
  
  keep <- names(outliers)[!outliers]
  pathway_metabolic_tpm <- pathway_metabolic_tpm[keep,]
  pathway_number_weight = 1 / gene_pathway_number[keep,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm, 1, function(x)by(x, all_cell_types, mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2, function(x) weighted.mean(x, pathway_number_weight/sum(pathway_number_weight)))
  mean_expression_shuffle[p, ] <-  mean_exp_pathway[cell_types]
  mean_expression_noshuffle[p, ] <-  mean_exp_pathway[cell_types]
    
  ##shuffle 5000 times:  
  ##define the functions 
  group_mean <- function(x){
    sapply(cell_types,function(y) rowMeans(pathway_metabolic_tpm[,shuffle_cell_types_list[[x]]==y,drop=F]))
  }
  column_weigth_mean <- function(x){
    apply(ratio_exp_eachCellType_list[[x]],2, function(y) weighted.mean(y, weight_values))
  }
  #####  
  times <- 1:5000
  weight_values <- pathway_number_weight/sum(pathway_number_weight)
  shuffle_cell_types_list <- lapply(times,function(x) sample(all_cell_types)) 
  names(shuffle_cell_types_list) <- times
  mean_exp_eachCellType_list <- lapply(times,function(x) group_mean(x))
  ratio_exp_eachCellType_list <- lapply(times,function(x) mean_exp_eachCellType_list[[x]] / rowMeans(mean_exp_eachCellType_list[[x]]))
  mean_exp_pathway_list <- lapply(times,function(x) column_weigth_mean(x))
  
  shuffle_results <- matrix(unlist(mean_exp_pathway_list),ncol=length(cell_types),byrow = T) 
  rownames(shuffle_results) <- times
  colnames(shuffle_results) <- cell_types
  for(c in cell_types){
    if(is.na(mean_expression_shuffle[p,c])) next
    if(mean_expression_shuffle[p,c]>1){
      pval <- sum(shuffle_results[,c] > mean_expression_shuffle[p,c]) / 5000 
    }else if(mean_expression_shuffle[p,c]<1){
      pval <- sum(shuffle_results[,c] < mean_expression_shuffle[p,c]) / 5000
    }
    if(pval>0.01) mean_expression_shuffle[p, c] <- NA  ### NA is  blank in heatmap
    pvalues_mat[p,c] <- pval
  }
}
all_NA <- rowAlls(is.na(mean_expression_shuffle))
mean_expression_shuffle <- mean_expression_shuffle[!all_NA,]
#heatmap
dat <- mean_expression_shuffle

sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1
# 2020-12-28 lly annoatote 
# lly think do not need to plot heatmap 
# pdf("/public/workspace/lily/res_cor_inte/Metabolic/KEGGpathway_activity_heatmap.pdf",onefile=T,width=6,height=9)
# mybreaks <- c(
#   seq(0.8, 0.85, length.out=23),
#   seq(0.86, 1.05, length.out=43),
#   seq(1.06, max(dat),length.out=34)
# ) 
# color <- colorRampPalette(c("blue","white","red"))(100)
# pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,
# color=color,breaks=mybreaks,fontsize_row=6,fontsize_col=6,cellwidth=9,cellheight=6)

# # default parmater
# mybreaks <- c(
#   seq(0, 0.5, length.out=33),
#   seq(0.51, 1.5, length.out=33),
#   seq(1.51, max(dat),length.out=34)
# ) 
# pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,
# color=color,breaks=mybreaks,fontsize_row=6,fontsize_col=6,cellwidth=9,cellheight=6)
# dev.off()

write.table(mean_expression_noshuffle,file="/public/workspace/lily/Lung2Brain/Version5/Metabolism/KEGGpathway_activity_noshuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(mean_expression_shuffle,file="/public/workspace/lily/Lung2Brain/Version5/Metabolism/KEGGpathway_activity_shuffle.txt",row.names=T,col.names=T,quote=F,sep="\t")
write.table(pvalues_mat,file="/public/workspace/lily/Lung2Brain/Version5/Metabolism/KEGGpathway_activity_shuffle_pvalue.txt",row.names=T,col.names=T,quote=F,sep="\t")

#boxplot show the distribution of pathway activity
# scRNA_dat <- as.data.frame(mean_expression_noshuffle)
# scRNA_dat$X <- NULL

# scRNA_df <- melt(scRNA_dat)
# scRNA_df <- scRNA_df[!is.na(scRNA_df$value),]
# g <- ggplot(scRNA_df,aes(x=variable,y=value,fill=variable)) +
#   scale_y_continuous(limits=c(0,3),breaks=0:3,labels=0:3)+
#   geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
#   stat_summary(fun.y = median,geom="point",size=1,color="blue")+
#   scale_fill_brewer(palette="Set2") +
#   theme_classic() + 
#   theme(legend.position="none",
#   axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
#   axis.text.y=element_text(colour="black", size = 6),
#   axis.line=element_line(size=0.2,color="black"),
#   axis.ticks = element_line(colour = "black",size=0.2),
#   panel.border = element_blank(), panel.background = element_blank(),
#   axis.ticks.length= unit(.5, "mm"))

# ggsave("/public/workspace/lily/res_cor_inte/Metabolic/pathway_activity_violinplot.pdf",g,width = 2.5,height=1.5,units="in",device="pdf",useDingbats=FALSE)










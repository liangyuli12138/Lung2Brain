#!/usr/bin/Rscript 
options(stringsAsFactors=F)

#===============================================================================================
# analysis TFs 
#===============================================================================================
# pyscenic 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tumor.dat <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type=="maliganant"))
# add module score to calculate BMS score 
#================================================================================================
tf.dat <- read.delim2("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant_auc_mtx.tsv",sep="\t")
rownames(tf.dat) <- tf.dat$Cell
tf.dat <- tf.dat[colnames(tumor.dat),]
tf.dat$Cell <- NULL
tf.dat <- apply(tf.dat,2,function(x){as.numeric(as.vector(x))})
rownames(tf.dat) <- colnames(tumor.dat)
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
tf.dat <- data.frame(tf.dat)
#================================================================================================
# use gene index to calculate result 
library(ineq)
gini_index <- apply(tf.dat, 2, function(x) {ineq(x, type = "Gini")})
tf.filter <- gini_index[which(gini_index[]>0.2)]
tf.dat <- tf.dat[,which(colnames(tf.dat)%in%names(tf.filter))]

# # calculate BMS score 
# source('~/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(tumor.dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
# mod <- as.data.frame(mod)
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_tumor_BMSscore.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_tumor_BMSscore.RDS")
# check result 
tf.dat$BMS <- mod$BMS_test_norm
res.tmp <- apply(tf.dat,2,function(x){
    tmp <- cor.test(as.numeric(as.vector(x)),as.numeric(as.vector(tf.dat$BMS)),method="spearman")
    percent <- length(which(x[]>0))/nrow(tf.dat)
    c(tmp$estimate,tmp$p.value,percent)
})

res.tmp <- t(res.tmp)
colnames(res.tmp) <- c("Cor","P.value","Percent")
res <- data.frame(res.tmp)
res$p.adj <- p.adjust(res$P.value)
res$log2p <- -log2(res$p.adj)

#=======================================================================
# 
# use TCGA RABIT result to intersect 
# Rabbit result 
#=======================================================================
load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData') # luad_mod
rabit_luad <- read.table('~/metastasis/data/TCGA_RABIT/LUAD/RABIT_LUAD.HiSeq.V2',sep='\t',header=T) #80 TFs 
rownames(rabit_luad) <- gsub("-",".",rabit_luad[,1])
rabit_luad$X <- NULL		#prepare

# get co-sample 
co_luad <- merge(luad_mod,rabit_luad,by='row.names')
rownames(co_luad) <- co_luad$Row.names
co_luad$Row.names <- NULL

# 2021-4-16 
# lly think maybe just use III and IV TCGA samaples, however is not good result 
# tmp <- read.delim("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_clinicalMatrix.txt",sep="\t",header=T)
# tcga.clin <- tmp[,c("sampleID","pathologic_M","pathologic_N","pathologic_T","pathologic_stage")]
# rownames(tcga.clin) <- gsub("-",".",tcga.clin$sampleID)
# co_luad.f <- co_luad[which(rownames(co_luad)%in%rownames(tcga.clin)[which(tcga.clin$pathologic_stage%in%c("Stage IIIA","Stage IIIB","Stage IV"))]),]


tf_luad <- data.frame(t(apply(co_luad[,-c(1:3)],2,function(x,y){tmp=cor.test(x,y);c(tmp$estimate,tmp$p.value,length(which(x!=0))/nrow(co_luad))},y=co_luad$BMS_test_norm)))
# col 1:3 is mod-result 
colnames(tf_luad) <- c('cor','pvalue','percent')
tf_luad$log2Pvalue <- -log2(tf_luad$pvalue) # -log2(pvalue)
tf_luad.f <- tf_luad[which(tf_luad$pvalue<0.05),]

#===========================================
gsub("\\.\\.\\.$","",rownames(res)) -> rownames(res)
res.f <- merge(tf_luad.f,res,by="row.names")

#===========================================
# try to find co_TFs
#===========================================
final.res <- res.f[which(res.f$cor*res.f$Cor>0),]
final.res <- final.res[which(final.res$p.adj<0.05),]
final.res$group <- "up"
final.res$group[which(final.res$cor<0)] <- "down"

#plot volcano plot 
#===========================================
library(ggplot2)
library(ggrepel)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TFs_correlation.pdf",useDingbats=F)
ggplot(final.res, aes(x = cor, y = log2Pvalue, colour=group,size=percent)) +
    geom_point(alpha=0.4) +
    scale_color_manual(values=c("#546de5","#ff4757"))+
    labs(x="Correaltion", y="-log2 (p-value)")+
    theme_bw()+
    geom_text(aes(y = log2Pvalue + .2, label = Row.names))
dev.off() 

saveRDS(final.res,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")






#===========================================
# plot 3D 


library(plot3D)
cex=1
phi=5
theta=10
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TFs_3D.pdf",useDingbats=F)

tmp <- (final.res[,c("cor","log2Pvalue","percent","Row.names")])    
scatter3D(tmp$cor, tmp$log2Pvalue, tmp$percent,pch=16,cex=cex,
    phi= phi,theta= theta,col=rgb(0.1,0.4,1,0.5))
#tmp <- (final.res[which(final.res$group=="down"),c("cor","log2Pvalue","percent")])    
# scatter3D(tmp$cor, tmp$log2Pvalue, tmp$percent,pch=16,cex=cex,
#     phi= phi,theta= theta,col="red",add=T)

text3D(tmp$cor, tmp$log2Pvalue, tmp$percent,labels=tmp$Row.names,add=T)

dev.off()

# plotly
library(plotly)
plot_ly(final.res, x = ~cor, y = ~log2Pvalue, z = ~percent, color = ~group , text = ~Row.names,size=1)










#=====================================================================================================================
# find traget genes 
#=====================================================================================================================
tf.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
up.tf <- tf.dat$Row.names[which(tf.dat$cor>0)]
tmp.res <- trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),]
up.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),"traget"]))
write.table(up.gene,file="./tmp_up_gene.txt",quote=F,row.names=F,col.names=F)
dn.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Repression"),"traget"]))
write.table(dn.gene,file="./tmp_dn_gene.txt",quote=F,row.names=F,col.names=F)










#==================================================================================================================
# transfactor TF MYBL2
#==================================================================================================================
file <- dir("~/metastasis/data/verify/GSE143145")[grep("\\.txt",dir("~/metastasis/data/verify/GSE143145"))]
dat.res <- matrix(NA)
for(i in 1:length(file)){
    tmp <- read.table(paste0("~/metastasis/data/verify/GSE143145/",file[i]),sep="\t",header=F)
    samplename <- strsplit(file[i],"_")[[1]][2]
    samplename <- gsub("-","_",samplename)
    colnames(tmp) <- c("ensembl_id",samplename)
    if(i==1){
        dat.res <- tmp 
    }else{
        dat.res <- merge(dat.res,tmp,by="ensembl_id")
    }
}

dat.res <- dat.res[-c(1:5),]
dat.res$ensembl_id <- gsub("\\..*","",dat.res$ensembl_id)
annol <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",sep="\t",header=T)
colnames(annol) <- c("ensembl_id","gene","gene.strat","gene.end")
final.res <- merge(dat.res,annol,by="ensembl_id")
final.res$genelen <- final.res$gene.end - final.res$gene.strat
tmp.res <- t(apply(final.res[,c(2:9,13)],1,function(x){
    tmp <- (x/x[9])*10^3
}))
res.f <- apply(tmp.res,2,function(x){
    tmp<- (x/sum(x))*10^6
})
res.f <- data.frame(res.f)
res.f$gene <- final.res$gene
res.f$genelen <- NULL 
aggregate(.~gene,data=res.f,FUN=median) -> res.ff
rownames(res.ff) <- res.ff$gene
res.ff$gene <- NULL

saveRDS(res.ff,file="/public/workspace/lily/metastasis/data/verify/GSE143145/GSE143145_NC_MYBL2_exp.RDS")

# #========================
# final.res$ensembl_id <- NULL
# final.res.f <- aggregate(.~gene,data=final.res,FUN=mean)

# rownames(final.res.f) <- final.res.f$gene
# final.res.f$gene <- NULL
# saveRDS(final.res.f,file="/public/workspace/lily/metastasis/data/verify/GSE143145/GSE143145_NC_MYBL2.RDS")

#====================================================================================================================================================
# prepare into expression file 
#====================================================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE143145/GSE143145_NC_MYBL2_exp.RDS")
# calculate ssgsea 
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=1000)
mod <- as.data.frame(mod)
mod$type <- c(rep("NC",4),rep("MYBL2",4))

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/GSE143145_BMS_mod.RDS")

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/MYBL2_kd_BMS.pdf")
boxplot(BMS_test_norm~type,data=mod)
legend('topright',legend=paste('p =',round(wilcox.test(BMS_test_norm~type,data=mod)$p.value,5)),bty='n')
dev.off()




#========================================================================================================================================
# sangji plot 
#========================================================================================================================================
# tumor2mac = read.table('results/malignant.Macrophage-pvalues.txt',header=T,row.names=1,sep='\t')
# tumor2mac.bin = apply(tumor2mac,2,function(x){ifelse(x<0.05,1,0)})
# colnames(tumor2mac.bin) = gsub('.malignant.Macrophage','',colnames(tumor2mac.bin))
# tumor2mac_lr = read.csv('results/nichenetr_result/malignant2Macrophage_WT_ligand_receptor.RDS',row.names = 1)
# tumor2mac_lt = read.csv('results/nichenetr_result/malignant2Macrophage_WT_ligand_target.RDS',row.names = 1)
# tumor2mac_lr.f = tumor2mac_lr#[which(tumor2mac_lr$weight>0.5),]
# tumor2mac_lt.f = tumor2mac_lt#[which(tumor2mac_lt$weight>0.004),]
# tumor2mac_lrt = merge(by='ligand',tumor2mac_lr,tumor2mac_lt.f,all.y=T)
# tumor2mac_lrt$flag = paste(tumor2mac_lrt$ligand,tumor2mac_lrt$receptor,sep='_')
# cm.lr = intersect(rownames(tumor2mac),tumor2mac_lrt$flag)
# tumor2mac_lrt.f = tumor2mac_lrt[which(tumor2mac_lrt$flag %in% cm.lr),]
# require(igraph)
# require(rCharts)
# require(rjson)
# require(plyr)
# require(reshape2)
# require(networkD3)
# nodes = data.frame(name=unique(c(unique(as.vector(tumor2mac_lrt.f$ligand)),
#                                  unique(as.vector(tumor2mac_lrt.f$receptor)),
#                                  unique(as.vector(tumor2mac_lrt.f$target)))))
# nodes$ID = seq(0,nrow(nodes)-1)
# link1 = unique(tumor2mac_lrt.f[,1:3])
# link2 = unique(tumor2mac_lrt.f[,c(2,4,5)])
# colnames(link1) = c("source", "target", "value")
# colnames(link2) = c("source", "target", "value")
# link1$value=link1$value*100
# link2$value=link2$value*100
# links = as.data.frame(rbind(link1,link2))
# links.map = merge(by.x='source',by.y='name',links,nodes)
# links.map = merge(by.x='target',by.y='name',links.map,nodes)
# sankeyNetwork(Links = links.map, Nodes = nodes,
#               Source = "ID.x", Target = "ID.y",
#               Value = "value", NodeID = "name",
#               fontSize= 12, nodeWidth = 30)

#=========================================================================================================================================
require(igraph)
require(rCharts)
require(rjson)
require(plyr)
require(reshape2)
require(networkD3)
options(stringsAsFactors=F)
#========================================================================================================================================
tf.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
up.tf <- tf.dat$Row.names[which(tf.dat$cor>0)]
tmp.res <- trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),]

#===========================================
# links 
links <- tmp.res[,c("source","traget","effect")]
links$effect <- 1
node <- unique(c(as.vector(tmp.res$source),as.vector(tmp.res$traget)))
node <- as.data.frame(node)

save(links,node,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/sangjiplot.RData")
 sankeyNetwork(Links = links, Nodes = node,
              Source = "source", Target = "traget",
              Value = "effect", NodeID = "node",
              fontSize= 12, nodeWidth = 30)



#=======================================================================================================================================
# use ggplot to plot sangji 
# 2020-11-27
# 2021-4-19
#=======================================================================================================================================
library(ggalluvial)
library(ggplot2)
library(RColorBrewer)
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/sangjiplot.RData")
tmp <- to_lodes_form(links,axes = 1:2,id = "Cohort") 

# plot 
cols <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75')
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Sangjiplot_TF.pdf",useDingbats=F,height=10)
ggplot(tmp,aes(x =factor(x,level = c("source","traget")),y=effect,stratum = stratum, alluvium = Cohort,fill = stratum, label =stratum)) +
geom_flow( width = 1/3) + # flow
geom_stratum( width = 1/3,linetype=0,size=0.5,alpha =0.5,color = "black")+ 
geom_text(stat ="stratum" , size =3) + #添加名字
scale_x_discrete(limits = c() )+ #去掉横坐标轴
theme_bw() +
theme(legend.position="none") +
scale_fill_manual(values = c(cols[1:4],rep("grey",110)))
dev.off()



# up gene KEGG pathway analysis 
#================================================================================
tf.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
up.tf <- tf.dat$Row.names[which(tf.dat$cor>0)]
tmp.res <- trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),]
up.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),"traget"]))


#######################################################################################################################
# 2021-4-26
# try to use cytoscape to run TFs and target gene
#======================================================================================================================
tf.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
trrust.dat.f <- trrust.dat[-which(trrust.dat$effect=="Unknown"),] # filter some unknow interaction

#==================================================
up.tf <- tf.dat$Row.names[which(tf.dat$cor>0)]
tmp.res <- trrust.dat.f[which(trrust.dat.f$source%in%up.tf),]
write.table(tmp.res,file="./up_TF.txt",row.names=F,col.names=T,quote=F,sep="\t")



















library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig

geneinfo <- bitr(up.gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.01)
tmp1 <- ekk@result
tmp1 <- tmp1[which(tmp1$pvalue<0.05),]





#==========================================================================================
# check if CEBPB and MYBL2 is differenta stage in fetal lung 
# 2020-12-10
# there is no CEBPB or MYBL2 in this data  
#==========================================================================================
# dat <- read.table("~/metastasis/data/verify/GSE121238/GSE121238_07Jun17_DA_filtered_annotated_RPKM.txt",sep="\t",header=T)
# tmp.f <- dat[,c(4,8:16,27)]
# # calculate TPM
# tmp <- t(apply(tmp.f[,2:11],1,function(x){x/x[length(x)]}))*10^3
# tmp.res <- apply(tmp,2,function(x){x/sum(x)})*10^6
# tmp.res <- data.frame(tmp.res)
# tmp.res$gene_name <- res$mouse_gene_name
# tmp.res$gene <- tmp.f$GeneName
# tmp.res$length <- NULL
# res.final <- aggregate(.~gene,data=tmp.res,FUN=median)
# rownames(res.final) <- res.final$gene
# res.final$gene <- NULL
# colnames(res.final) <- c(paste0("10.3w_",1:3),paste0("14.5w_",1:2),paste0("14.4w_",1),paste0("20.4w_",1),paste0("20.5w_",1:2))
# saveRDS(res.final,file="/public/workspace/lily/metastasis/data/verify/GSE121238/GSE121238_exp.RDS")



#========================================================================================
# other data 
# GSE14334 data
#========================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14334/GSE14334_series_matrix.txt",sep="\t",header=T,comment.char="!")
gpl570 <- read.table("/public/workspace/lily/metastasis/data/verify/GPL570.txt",sep="\t",header=T)
# merge result 
tmp.res <- merge(gpl570,dat,by.x="ID",by.y="ID_REF")
tmp.res.f <- tmp.res[-grep("///",tmp.res$Gene.Symbol),]
tmp.res.f$ID <- NULL
tmp.res.f$Gene.Symbol <- as.vector(tmp.res.f$Gene.Symbol)
aggregate(.~Gene.Symbol,data=tmp.res.f,FUN=median) -> res.final
rownames(res.final) <- as.vector(res.final$Gene.Symbol)
res.final$Gene.Symbol <- NULL
# do some filter ,1:25 is unknow gene 
res.final.f <- res.final[-(1:25),]
ann <- read.table("~/metastasis/data/verify/GSE14334/GSE14334_ann.txt.txt",sep="\t")
colnames(res.final.f) <- ann$V2
a <- res.final.f[,order(colnames(res.final.f))]
a <- a[,c(21:38,1:20)]
saveRDS(a,file="/public/workspace/lily/metastasis/data/verify/GSE14334/GSE1433_expr.RDS")

#=======================================================================================
# use denssity plot to show result 
# 2020-12-11
#=======================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14334/GSE1433_expr.RDS")
res <- as.data.frame(t(dat[c("MYBL2","CEBPB"),]))
res$type <- sapply(strsplit(rownames(res),"_"),function(x){paste0(x[1],x[2])})
tmp.res <- aggregate(.~type,data=res,FUN=median)
tmp.f <- melt(tmp.res)
colnames(tmp.f)[2] <- c("Gene")
tmp.f$type <- factor(as.vector(tmp.f$type),levels=c(unique(tmp.f$type)[c(14:29,1:13)]))
tmp.f$rough.type <- "unknow"
tmp.f$rough.type[which(tmp.f$type%in%c("lung53","lung59","lung63"))] <- "Day_50-60"
tmp.f$rough.type[which(tmp.f$type%in%c("lung72","lung74","lung76","lung78"))] <- "Day_70-80"
tmp.f$rough.type[which(tmp.f$type%in%c("lung82","lung83","lung85","lung87","lung89"))] <- "Day_80-90"
tmp.f$rough.type[which(tmp.f$type%in%c("lung91","lung94","lung96","lung98"))] <- "Day_90-100"
tmp.f$rough.type[which(tmp.f$type%in%c("lung101","lung103","lung105","lung108"))] <- "Day_100-110"
tmp.f$rough.type[which(tmp.f$type%in%c("lung110","lung113","lung117"))] <- "Day_110-120"
tmp.f$rough.type[which(tmp.f$type%in%c("lung122","lung130","lung133","lung134"))] <- "Day_122-140"
tmp.f$rough.type[which(tmp.f$type%in%c("lung140","lung154"))] <- "Day_140-154"

tmp.f$rough.type <- factor(tmp.f$rough.type,levels=c("Day_50-60","Day_70-80","Day_80-90","Day_90-100","Day_100-110","Day_110-120","Day_122-140","> y_140-154"))

# ggplot 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/Fetal_lung_MYBL2_CEBPB.pdf",useDingbats=F,width=12)
ggplot(data=tmp.f,aes(x=rough.type,y=value,fill=Gene))+geom_boxplot() + theme(panel.grid.major = element_blank()) +
theme_classic()+labs(x="Stage",y='Gene Expression value')
dev.off()



#=========================================================================================
# ggplot GSE131907 data 
# 2020-12-15
# use sangsung data to verify 
#=========================================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat <- subset(dat,cells=which(dat$Cell_subtype=="Malignant cells"&dat$Sample_Origin=="mBrain"))
sub.dat <- recluster(sub.dat)

#======================
# 
#======================
data <- sub.dat[['RNA']]@data
sub.dat$type.TF <- "none"
sub.dat$type.TF[which(data["CEBPB",]>0&data["MYBL2",]>0)] <- "Both"
sub.dat$type.TF[which(data["CEBPB",]==0&data["MYBL2",]>0)] <- "MYBL2"
sub.dat$type.TF[which(data["CEBPB",]>0&data["MYBL2",]==0)] <- "CEBPB"

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE131907_MYBL2_CEBPB.pdf",useDingbats=F)
DimPlot(sub.dat,group.by="type.TF",cols=c("#ed1b2e","#ecb731","#0091cd","#d7d7d8"))
#=======================
num.both <- length(which(data["CEBPB",]>0&data["MYBL2",]>0))/ncol(data)
num.CEBPB <- length(which(data["CEBPB",]>0&data["MYBL2",]==0))/ncol(data)
num.MYBL2 <- length(which(data["CEBPB",]==0&data["MYBL2",]>0))/ncol(data)
barplot(c(num.both,num.CEBPB,num.MYBL2),ylim=c(0,0.5),names=c("Both","CEBPB","MYBL2"),col=c("#ed1b2e","#ecb731","#0091cd"))
dev.off()


#==========================================================================
# difference in LC and Brain metastasis 
# 2020-12-18
#==========================================================================
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE131907_LC_BM_MYBL2_CEBPB_percent.pdf")
barplot(c(1019/15423,314/7270),ylim=c(0,0.08),names=c("BM","LC"))
legend("topright",legend=paste0("P= ",fisher.test(matrix(c(314,1019,7270,15423),nrow=2))$p.value))
dev.off()




#=============================
# make sure early Lung is no both 
#=============================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat <- subset(dat,cells=which(dat$Cell_type=="Epithelial cells"&dat$Sample_Origin=="tLung"))
#sub.dat <- recluster(sub.dat)
data <- sub.dat[['RNA']]@data

# another try 
# 2020-12-16
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat <- subset(dat,cells=which(dat$Cell_type=="Epithelial cells"&dat$Sample_Origin=="tL/B"))
#sub.dat <- recluster(sub.dat)
data <- sub.dat[['RNA']]@data
num.both <- length(which(data["CEBPB",]>0&data["MYBL2",]>0))/ncol(data)
num.CEBPB <- length(which(data["CEBPB",]>0&data["MYBL2",]==0))/ncol(data)
num.MYBL2 <- length(which(data["CEBPB",]==0&data["MYBL2",]>0))/ncol(data)




#==========================================================================================
# check another lung adenocarcinoma data 
# 2020-12-16
# GSE117570 
#==========================================================================================
library(Seurat)
file <- dir("/public/workspace/lily/metastasis/data/verify/GSE117570")
for(i in 1:length(file)){
    data <- read.table(paste0("/public/workspace/lily/metastasis/data/verify/GSE117570/",file[i]),sep="\t",header=T)
    num.both <- length(which(data["CEBPB",]>0&data["MYBL2",]>0))/ncol(data)
    num.CEBPB <- length(which(data["CEBPB",]>0&data["MYBL2",]==0))/ncol(data)
    num.MYBL2 <- length(which(data["CEBPB",]==0&data["MYBL2",]>0))/ncol(data)
    print(c(num.both,num.CEBPB,num.MYBL2))
}






#========================================================================
# CEBPB and MYBL2 bin to show 
# 2020-12-17
# change colour 
#========================================================================
library(monocle)
library(Seurat)

cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/MYBL2_CEBPB_tumor_trajectory.pdf",useDingbats=F)
plot_cell_trajectory(cds,color_by="type_group")+scale_colour_manual(values=c("#205527","#0863b5"))
cds$type.TF <- "unknow"
cds$type.TF[which(cds@assayData$exprs["CEBPB",]>0&cds@assayData$exprs["MYBL2",]==0)] <- "CEBPB"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]==0)] <- "MYBL2"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]>0)] <- "Both"
plot_cell_trajectory(cds,color_by="type.TF")+scale_colour_manual(values=c("#ffbd0a","#46b7fd","#7ac142"))
dev.off()


# calculate percentage
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
cds$type.TF[which(cds@assayData$exprs["CEBPB",]>0&cds@assayData$exprs["MYBL2",]==0)] <- "CEBPB"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]==0)] <- "MYBL2"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]>0)] <- "Both"

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_tumor_MYBL2_CEBPB_percentage.pdf")
#table(cds$type.TF,cds$type_group)
#table(dat$type,dat$type_group)
barplot(c(28/2078,908/8243),names=c("LC","BM"),ylim=c(0,0.15))
dev.off()


# calculate if Hallmarks is different in 3 groups
#=========================================================================
cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
data <- as.matrix(cds@assayData$exprs)
source('~/software/ssGSEA/ssgseaMOD.r')
modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
mod <- mod.analyze2(data,modlist,'/public/workspace/lily/MOD_file/HALLMARK/',permN=0)
mod <- data.frame(mod)
cds$type.TF <- "unknow"
cds$type.TF[which(cds@assayData$exprs["CEBPB",]>0&cds@assayData$exprs["MYBL2",]==0)] <- "CEBPB"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]==0)] <- "MYBL2"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]>0)] <- "Both"
mod$type <- cds$type.TF
mod.f <- mod[,c(51:100,151)]
saveRDS(mod.f,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_LC_TF_Hallmark_mod.RDS")
mod.f <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/LCBM_LC_TF_Hallmark_mod.RDS")
tmp <- data.frame(aggregate(.~type,data=mod.f,FUN=mean))
rownames(tmp) <- tmp$type
tmp$type <- NULL
tmp.res <- t(tmp)
rownames(tmp.res) <- gsub("^HALLMARK_|_norm$","",rownames(tmp.res))
rownames(tmp.res) <- gsub("_"," ",rownames(tmp.res))

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_LC_TF_hallmark.pdf")
pheatmap::pheatmap(tmp.res,scale="row")
dev.off()




# 2021-1-30 
# re-plot heatmap

plot_Circ_heatmap <- function(mat, cluster) {
    # ref: https://zhuanlan.zhihu.com/p/136138642
    #@mat: row(sample or groups) X col(pair info or genes)
    #@cluster: set cluster number in hclust function
    library(dendextend)
    library("circlize")
    library(RColorBrewer)

    mat=scale(mat, center = TRUE, scale = TRUE)
    dend <-as.dendrogram(hclust(dist(t(mat))))
    n=1
    dend <-dend %>% set("branches_k_color", k = n) 
    #par(mar=c(7.5,3,1,0))
    #plot(dend)
    # 聚类后的样本信息
    mat2 = mat[, order.dendrogram(dend)]
    lable1=row.names(mat2)
    lable2=colnames(mat2)
    nr = nrow(mat2)
    nc = ncol(mat2)
    col_fun = colorRamp2(c(-1.5, 0, 1.5), c("skyblue", "white", "red"))
    col_mat = col_fun(mat2)
    par(mar=c(0,0,0,0))
    circos.clear()
    circos.par(
        canvas.xlim =c(-2,2),
        canvas.ylim = c(-2,2),
        cell.padding = c(0,0,0,0), 
        gap.degree =90
    )
    factors = "a"
    circos.initialize(factors, xlim = c(0, ncol(mat2)))
    circos.track(
        ylim = c(0, nr),bg.border = NA,track.height = 0.1*nr,
        panel.fun = function(x, y) {
            for(i in 1:nr) {
                circos.rect(xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
                    xright = 1:nc, ytop = rep(nr - i + 1, nc),
                    border = "black",
                    col = col_mat[i,]
                )
                circos.text(x = nc,
                    y = 6.4 -i,
                    labels = lable1[i],
                    facing = "downward", niceFacing = TRUE,
                    cex = 0.6,
                    adj = c(-0.2, 0))
            }
        }
    )
    for(i in 1:nc){
        circos.text(x = i-0.4,
        y = 5,
        labels = lable2[i],
        facing = "clockwise", niceFacing = TRUE,
        cex = 0.4,adj = c(0, 0))
    }
    #添加树
    max_height <-max(attr(dend, "height"))
    circos.track(ylim = c(0, max_height),bg.border = NA,track.height = 0.3, 
        panel.fun = function(x, y){
        circos.dendrogram(dend = dend,
        max_height = max_height)
    })
    circos.clear()
    # 添加图例
    library(ComplexHeatmap)
    lgd <- Legend(at = c(-2,-1, 0, 1, 2), col_fun = col_fun, 
                title_position = "topcenter",title = "Z-score")
    draw(lgd, x = unit(0.7, "npc"), y = unit(0.7, "npc"))
}

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LCBM_LC_TF_hallmark_circle.pdf")
plot_Circ_heatmap(t(tmp.res))
dev.off()




#=================================================================================
# calculate metabolism in 3groups 
# 2020-12-26
#=================================================================================
cds <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
cds$type.TF <- "unknow"
cds$type.TF[which(cds@assayData$exprs["CEBPB",]>0&cds@assayData$exprs["MYBL2",]==0)] <- "CEBPB"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]==0)] <- "MYBL2"
cds$type.TF[which(cds@assayData$exprs["MYBL2",]>0&cds@assayData$exprs["CEBPB",]>0)] <- "Both"

saveRDS(cds,file="/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")

data <- as.matrix(cds@assayData$exprs)
source('~/software/ssGSEA/ssgseaMOD.r')




















#========================================================================================
# use GSE143423 data to verify 
# 2020-12-15
#========================================================================================
#====================================================================================
# use GSE143423 to verfiy CD4 and ANXA1 T cell 
# FeaturePlot
#====================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE143423/GSE143423_seurat.RDS")
# define cell type 
#========================================
pdf("./tmp.pdf")
FeaturePlot(dat,features=c('CD3D','CD3E','CD2','PTPRC'),label=T,order=F) # T cell
FeaturePlot(dat,features=c('CD19','CD68','FCGR3A','FCGR1A'),label=T,order=F) # Myeloid
FeaturePlot(dat,features=c('MS4A1',"CD79A",'PTPRC'),label=T,order=T) # B cell 
FeaturePlot(dat,features=c('MAG','MOG','CNDP1','PTPRC'),label=T,order=T) # Oligodendrocyte
FeaturePlot(dat,features=c('COL1A1','COL1A2','DCN','CD248'),label=T,order=T) # Fibroblast/Vascular
FeaturePlot(dat,features=c('CLDN5','VWF','ABCG2','CDH5'),label=T,order=T) # Endothelial
FeaturePlot(dat,features=c("EGFR","EPCAM","PTPRC"),label=T,order=T,cols=c("lightgrey", "red"))
dev.off()
tumor <- subset(dat,cells=which(dat$seurat_clusters%in%c(0,12,3,4,5,7,9)))
data <- tumor[['RNA']]@data
num.both <- length(which(data["CEBPB",]>0&data["MYBL2",]>0))/ncol(data)
num.CEBPB <- length(which(data["CEBPB",]>0&data["MYBL2",]==0))/ncol(data)
num.MYBL2 <- length(which(data["CEBPB",]==0&data["MYBL2",]>0))/ncol(data)
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE143423_MYBL2_CEBPB.pdf",useDingbats=F)
barplot(c(num.both,num.CEBPB,num.MYBL2),ylim=c(0,0.5),names=c("Both","CEBPB","MYBL2"),col=c("#ed1b2e","#ecb731","#0091cd"))
dev.off()











#==============================================================================================================
# 2020-12-28
# check MYBL2 and CEBPB target genes 
#==============================================================================================================
tf.dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
up.tf <- tf.dat$Row.names[which(tf.dat$cor>0)]
tmp.res <- trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect!="Unknown"),]

# check CEBPB activation genes 
cebpb_ac_gene <-  as.vector(tmp.res[which(tmp.res$source=="CEBPB"&tmp.res$effect=="Activation"),"traget"])
cebpb_rp_gene <-  as.vector(tmp.res[which(tmp.res$source=="CEBPB"&tmp.res$effect=="Repression"),"traget"])

# KEGG pathway 
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig
geneinfo <- bitr(cebpb_rp_gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.01)
tmp1 <- ekk@result
cebpb_rp_pathway <- tmp1[which(tmp1$pvalue<0.05),]



# check MYBL2 activation genes 
mybl2_ac_gene <-  as.vector(tmp.res[which(tmp.res$source=="MYBL2"&tmp.res$effect=="Activation"),"traget"])
mybl2_rp_gene <-  as.vector(tmp.res[which(tmp.res$source=="MYBL2"&tmp.res$effect=="Repression"),"traget"])

# KEGG pathway 
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig
geneinfo <- bitr(mybl2_ac_gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.01)
tmp1 <- ekk@result
mybl2_ac_pathway <- tmp1[which(tmp1$pvalue<0.05),]





up.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%up.tf&trrust.dat$effect=="Activation"),"traget"]))


















#======================================================================================================================
# Analysis metabolism 
# 2020-12-28
#======================================================================================================================
library(pheatmap)
dat <- read.table("/public/workspace/lily/Lung2Brain/inte7/metabolism/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
dat[is.na(dat)] <- 1
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/MYBL2_CEBPB_tumor_metabolism.pdf",height=15,width=10)
pheatmap::pheatmap(dat,scale="row")
dev.off()


# 2021-1-30
dat <- read.table("/public/workspace/lily/Lung2Brain/inte7/metabolism/KEGGpathway_activity_shuffle.txt",sep="\t",header=T)
dat[is.na(dat)] <- 1
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/MYBL2_CEBPB_tumor_metabolism_circle.pdf",height=15,width=10)
plot_Circ_heatmap(t(dat))
dev.off()









#====================================================================================================================
# check Both is HPSC 
# 2021-1-29 LLY check this data 
# 尝试通过使用文章中 HPSC 的signature 对双阳的细胞进行分析，但是没有找到准确的signature，可能可以使用C5 的差异基因进行计算
# 1. 使用C5的差异基因进行计算，找到人的同源基因，然后进行分析
# 2. 查看该数据在Lung cancer evolution的过程中TF的表达值的变化
#===================================================================================================================
# 1. check C5 siganture varation 
tmp.gene <- as.vector(read.table("/public/workspace/lily/metastasis/data/verify/LUAD_evolution_mm/C5_markergene.txt",sep="\t",header=F)[,1])
tmp.dat <- readRDS("/public/workspace/lily/metastasis/data/verify/Human_mm_gene.RDS")
genelist <- tmp.dat[which(tmp.dat$external_gene_name%in%tmp.gene),"hsapiens_homolog_associated_gene_name"]
source('~/software/ssGSEA/ssgseaMOD.r')
mod.generate(genelist,'HPSC_C5',out='/public/workspace/lily/MOD_file/HPSC_C5.mod')

# calculate mod
setwd("/public/workspace/lily/Lung2Brain/")
source('~/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/trajectory/LCBM_LC_tumor_trajectory_MYBL2_CEBPB.RDS")
mod <- mod.analyze2(as.matrix(dat@assayData$exprs),c('HPSC_C5'),'/public/workspace/lily/MOD_file/',permN=0)
mod <- as.data.frame(mod)

mod$type <- dat$type.TF

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/HPSC_MYBL2_CEBPB.pdf",useDingbats=F)
boxplot(HPSC_C5_norm~type,data=mod,FUN=median,outline=F,ylim=c(0,1))
dev.off()





# 2. use data to caculate 
#===================================================================================================================
library(hdf5r)

tmp <- H5File$new("./GSE154989_mmLungPlate_fQC_dSp_rawCount.h5",mode="r")

















































############################################################################################################################################################
############################################################################################################################################################
# 2021-2-22
# calculate MYBL2 CEBPB Both and None type tumor cell in scRNA LC samples 
# if this group have different BMS score 
#===========================================================================================================================================================
# 1. SanSung Group 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_tumor.RDS")

# set cell type 
sub.dat <- subset(dat,cells=which(dat$Sample_Origin=="tLung"))
sub.dat$celltype <- "unknow"
tmp.data <- sub.dat[["RNA"]]@data
sub.dat$celltype[tmp.data["CEBPB",]>0&tmp.data["MYBL2",]==0] <- "CEBPB"
sub.dat$celltype[tmp.data["CEBPB",]==0&tmp.data["MYBL2",]>0] <- "MYBL2"
sub.dat$celltype[tmp.data["CEBPB",]==0&tmp.data["MYBL2",]==0] <- "None"
sub.dat$celltype[tmp.data["CEBPB",]>0&tmp.data["MYBL2",]>0] <- "Both"



# load BMS gene 
# use addModulescore
gene <- readRDS("/public/workspace/lily/Lung2Brain/inte7/BMS_gene.RDS")
sub.dat <- AddModuleScore(sub.dat,features=list(gene))

# use ssgsea
# show good result 
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[["RNA"]]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- as.data.frame(mod)
mod$type <- factor(sub.dat$celltype,levels=c("Both","MYBL2","CEBPB","None"))
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/GSE131907_Lungtumor_TF_BMS_mod.RDS")
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/GSE131907_Lungtumor_TF_BMS.pdf")
boxplot(BMS_test_norm~type,data=mod,FUN=median,outline=F)
dev.off()





#=======================================================================================================================================================
# 2. calculate in self data 
##########################################################################################################
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$type=="maliganant"&dat$type_group=="LC"))

tmp.data <- sub.dat[['RNA']]@data
sub.dat$TFgroup <- "unknow"
sub.dat$TFgroup[tmp.data["CEBPB",]>0&tmp.data["MYBL2",]==0] <- "CEBPB"
sub.dat$TFgroup[tmp.data["CEBPB",]==0&tmp.data["MYBL2",]>0] <- "MYBL2"
sub.dat$TFgroup[tmp.data["CEBPB",]==0&tmp.data["MYBL2",]==0] <- "None"
sub.dat$TFgroup[tmp.data["CEBPB",]>0&tmp.data["MYBL2",]>0] <- "Both"

# calculate BMS score 
source('~/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(sub.dat[["RNA"]]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- as.data.frame(mod)
mod$TFgroup <- factor(sub.dat$TFgroup,levels=c("Both","MYBL2","CEBPB","None"))

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/LC_TF_BMS.pdf")
boxplot(BMS_test_norm~TFgroup,data=mod,FUN=median,outline=F,ylim=c(0,1))
dev.off()
















































##################################################################################################################################
# 2021-3-12
# Feature plot of all malignant cells with MYBL2 and CEBPB 
# TF activity
#=================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tumor <- subset(dat,cells=which(dat$type=="maliganant"))
# tf.dat <- read.delim2("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_auc_mtx.tsv",sep="\t")
# rownames(tf.dat) <- tf.dat$Cell
# tf.dat <- tf.dat[colnames(tumor),]

# tumor$MYBL2_ac <- as.numeric(as.vector(tf.dat[,"MYBL2..."]))
# tumor$CEBPB_ac <- as.numeric(as.vector(tf.dat[,"CEBPB..."]))
tmp.dat <- tumor[['RNA']]@data
tumor$MYBL2 <- tmp.dat["MYBL2",]
tumor$CEBPB <- tmp.dat["CEBPB",]

pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/CEBPB_cutoff_featureplot.pdf",useDingbats=F)
FeaturePlot(tumor,features="CEBPB",min.cutoff=1,order=T)
FeaturePlot(tumor,features="MYBL2",min.cutoff=0,order=T)
dev.off()









# 2021-4-15
# test in Multiple Lung brain sample 
# and in  ILCM sample 
#===================================================================================================================================
library(Seurat)
D0927 <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/D0927.RDS")
tumor.D <- subset(D0927,cells=which(D0927$infercnv.type=="Tumor"))
saveRDS(tumor.D,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
E0927 <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/E0927.RDS")
tumor.E <- subset(E0927,cells=which(E0927$infercnv.type=="Tumor"))
saveRDS(tumor.E,file="/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")

# run pyscenic in APP 
# /public/workspace/mirrors/scenic_ref/human
# and also use MSK data 
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE123904/GSE123902.RDS")
sub.dat <- subset(dat,cells=which(dat$seurat_clusters==6&dat$group=="PRIMARY"))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")


# and maybe should check 
# 2021-4-19
#####################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$maliganant=="tumor"&dat$type_group=="LC"))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/inte7_LC/LC.tumor.RDS")























#####################################################################################################################################
# 2021-4-19
# check result for MSK data 
# 0. GSE123902
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")
tf.dat <- read.csv("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/step3.auc_mtx.csv")
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))

# use BMS score 
# sub.dat <- AddModuleScore(sub.dat,features=gene,name="BMS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
mod$group <- "Early"
mod$group[which(dat$sample=="LX675")] <- "Advanced"
######################################################################################################################################
mod$MYBL2 <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
mod$CEBPB <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
mod$BACH1 <- as.numeric(as.vector(tf.dat[,"BACH1"]))

# check result is OK , however do not know how to show result 



##################################################################################################################################
# 1. Multiple lung cancer show result 
#=================================================================================================================================
# 1.1 D0927
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
tf.dat <- read.csv("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/D0927/step3.auc_mtx.csv")
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)
######################################################################################################################################
mod$MYBL2 <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
mod$CEBPB <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
mod$BACH1 <- as.numeric(as.vector(tf.dat[,"BACH1"]))





# 1.2 E0927
#====================================================================================================================================
#
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
tf.dat <- read.csv("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/E0927/step3.auc_mtx.csv")
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
tf.dat <- tf.dat[which(rownames(tf.dat)%in%colnames(dat)),]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)

mod$MYBL2 <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
mod$CEBPB <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
mod$BACH1 <- as.numeric(as.vector(tf.dat[,"BACH1"]))





######################################################################################################################################
# 2. inte7  Lung cancer 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/inte7_LC/LC.tumor.RDS")
tf.dat <- read.csv("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/inte7_LC/step3.auc_mtx.csv")
rownames(tf.dat) <- tf.dat$Cell
tf.dat$Cell <- NULL
colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
mod <- data.frame(mod)

mod$MYBL2 <- as.numeric(as.vector(tf.dat[,"MYBL2"]))
mod$CEBPB <- as.numeric(as.vector(tf.dat[,"CEBPB"]))
mod$BACH1 <- as.numeric(as.vector(tf.dat[,"BACH1"]))







# 3. use Bulk RNA seq to check TFs 
# 2021-4-20
# Run in R-4.0.2
#=================================================================================================================
library(dorothea)
library(bcellViper)
library(dplyr)
library(viper)
data(dorothea_hs,package="dorothea")
#subset DoRothEA to the confidence levels A and B to include only the high quality regulons
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE72094/GSE72094_expr.RDS")
regulons = dorothea_hs[which(dorothea_hs$confidence%in%c("A","B")),]
tf_activities <- run_viper(dat, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))

rs <- as.data.frame(t(tf_activities[c("MYBL2","CEBPB"),]))
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE72094/KRAS_info.txt",header=T))
ann <- ann[-1,,drop=F]
colnames(ann) <- "KRAS_type"
res.final <- merge(rs,ann,by="row.names")











# 2021-4-15
# 2021-4-20 update, lly think should use Brain metastasis sample to do this 
# calculate TFs activaity in all malignant cell and BMS score correlation 
#==================================================================================================================================
# dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
# tumor <- subset(dat,cells=which(dat$type=="maliganant"))
# #tumor <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type=="maliganant"))
# tf.dat <- read.delim2("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_auc_mtx.tsv",sep="\t")
# rownames(tf.dat) <- tf.dat$Cell
# tf.dat <- tf.dat[colnames(tumor),]
# tf.dat$Cell <- NULL
# tf.dat <- apply(tf.dat,2,function(x){as.numeric(as.vector(x))})
# rownames(tf.dat) <- colnames(tumor)
# colnames(tf.dat) <- gsub("\\.\\.\\.$","",colnames(tf.dat))
# tf.dat <- data.frame(tf.dat)
# # calculate BMS score use ssgsea
# #==================================================================================================================================
# source('~/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(tumor[['RNA']]@data),c('BMS_test'),'/public/workspace/lily/Lung2Brain/inte7/',permN=0)
# mod <- as.data.frame(mod)

# tf.dat$BMS <- mod$BMS_test_norm
# tf.dat$type_group <- tumor$type_group

# res.tmp <- apply(tf.dat,2,function(x){
#     tmp <- cor.test(as.numeric(as.vector(x)),as.numeric(as.vector(tf.dat$BMS)))
#     c(tmp$estimate,tmp$p.value)
# })

# res.tmp <- t(res.tmp)
# colnames(res.tmp) <- c("Cor","P.value")
# res <- data.frame(res.tmp)
# res$p.adj <- p.adjust(res$P.value)
# res$log2p <- -log2(res$p.adj)

# # use TCGA rabit result to filter 
# load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData') # luad_mod
# rabit_luad <- read.table('~/metastasis/data/TCGA_RABIT/LUAD/RABIT_LUAD.HiSeq.V2',sep='\t',header=T) #80 TFs 
# rownames(rabit_luad) <- gsub("-",".",rabit_luad[,1])
# rabit_luad$X <- NULL		#prepare

# # get co-sample 
# co_luad <- merge(luad_mod,rabit_luad,by='row.names')
# rownames(co_luad) <- co_luad$Row.names
# co_luad$Row.names <- NULL
# tf_luad <- data.frame(t(apply(co_luad[,-c(1:3)],2,function(x,y){tmp=cor.test(x,y);c(tmp$estimate,tmp$p.value,length(which(x!=0))/nrow(co_luad))},y=co_luad$BMS_test_norm)))
# # col 1:3 is mod-result 
# colnames(tf_luad) <- c('cor','pvalue','percent')
# tf_luad$log2Pvalue <- -log2(tf_luad$pvalue) # -log2(pvalue)
# tf_luad.f <- tf_luad[which(tf_luad$pvalue<0.05),]

# res.f <- merge(tf_luad.f,res,by="row.names")

# #===========================================
# #===========================================
# # try to find co_TFs
# #===========================================
# final.res <- res.f[which(res.f$cor*res.f$Cor>0),]
# final.res <- final.res[which(final.res$p.adj<0.05),]
# final.res$group <- "up"
# final.res$group[which(final.res$cor<0)] <- "down"

# #plot volcano plot 
# #===========================================
# library(ggplot2)
# library(ggrepel)
# pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TFs_correlation.pdf",useDingbats=F)
# ggplot(final.res, aes(x = cor, y = log2Pvalue, colour=group,size=percent)) +
#     geom_point(alpha=0.4) +
#     scale_color_manual(values=c("#546de5","#ff4757"))+
#     labs(x="Correaltion", y="-log2 (p-value)")+
#     theme_bw()+
#     geom_text(aes(y = log2Pvalue + .2, label = Row.names))
# dev.off() 

# saveRDS(final.res,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")


# use scatter plot 3d to show result 
# 2021-4-19
#====================================================================================================================
library(scatterplot3d)

final.res <- readRDS("/public/workspace/lily/Lung2Brain/inte7/res_Data/TFs_cor.RDS")
tmp <- data.frame(final.res[,c("cor","log2Pvalue","percent","Row.names")])    
pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/TF_scatter3d.pdf",useDingbat=F)
colors<-ifelse(tmp$cor<0,"blue","red")
scatterplot3d(x=tmp$cor,y=tmp$log2Pvalue,z=tmp$percent, pch = 16, angle=10,
              color =colors,xlab="Correlation",ylab="log2pvalue",zlab="Percentage",type="h"
    )

dev.off()







#########################################################################################################################################################
# 2021-4-20 
# run Pyscienc for 3 sample : D0927 E0927 and MSK  # have done in 2021-4-19 
# and now need to use GSE
#========================================================================================================================================================
library(Seurat)
# D0927 
dat <- readRDS()






























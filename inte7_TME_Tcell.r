
#===============================================================================
# Lung2Brain T cell  analysis 
# 
#===============================================================================
# 1. TME analysis should add GBM data 
# 2. T cell shoule classify into subset 
#===============================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
T_cell <- subset(dat,cells=which(dat$type=="T_cell"))


inte.list <- list() 
samplelist <- unique(T_cell$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(T_cell,cells=which(T_cell$orig.ident==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=20)
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
inte <- RunUMAP(inte,dims=1:10)

saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
#===============================================================================
# re clustering
#===============================================================================
# DefaultAssay(T_cell) <- "integrated"
# T_cell <- FindVariableFeatures(object = T_cell)
# # scaling
# all.genes <- rownames(x = T_cell)
# T_cell <- ScaleData(object = T_cell, features = all.genes)
# # PCA
# T_cell <- RunPCA(object = T_cell, features = VariableFeatures(object = T_cell))
# # clustering
# T_cell <- FindNeighbors(object = T_cell,dims=1:10)
# # select proper resolution
# T_cell <- FindClusters(object = T_cell)
# # T-SNE
# T_cell <- RunTSNE(object = T_cell,dims=1:10,check_duplicates = FALSE)
# T_cell <- RunUMAP(T_cell,dims=1:10)

#saveRDS(T_cell,file="/public/workspace/lily/Lung2Brain/TME/T_cell/inte_T_cell.RDS")




#===============================================================================
# run cellphoneDB in MDM cells and T cells
# do not get subset first 
# 2020-12-10 
#===============================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
Tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
# set refine cell type 
dat$type.CCC <- dat$type
dat$type.CCC[which(colnames(dat)%in%colnames(Myeloid))] <- Myeloid$type.refine
dat$type.CCC[which(colnames(dat)%in%colnames(Tcell))] <- Tcell$type.rough
# subset LCBM cells 
sub.dat <- subset(dat,cells=which(dat$type.CCC%in%c("BMDM","CD4_naive","CD4_Tcell","CD8_exhausted","CD8_Tcell","Treg")&dat$type_group=="LCBM"))
sub.dat$sample <- sub.dat$orig.ident
sub.dat$sample[grep("RD-20180817-001-SR18271",sub.dat$sample)] <- "lesion1"
sub.dat$sample[grep("RD-20180817-002-SR18271",sub.dat$sample)] <- "lesion2"
sub.dat$sample[grep("T-Bsc1",sub.dat$sample)] <- "T_Bsc1"


cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}

sample_name <- unique(sub.dat$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/",sample_name[i],"/")
    tmp <- subset(sub.dat,cells=which(sub.dat$sample==sample_name[i]))
    cellphoneDB_input(tmp,"type.CCC",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

#=======================================================================================================================
# run in shell 
for i in `ls /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done








#======================================================================================================================
# do cellphoneDB with Endothelials and MDMs with CD4 + Treg
# 2020-12-12
#======================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
Tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
# set refine cell type 
dat$type.CCC <- dat$type
dat$type.CCC[which(colnames(dat)%in%colnames(Myeloid))] <- Myeloid$type.refine
dat$type.CCC[which(colnames(dat)%in%colnames(Tcell))] <- Tcell$type.rough
# subset LCBM cells 
sub.dat <- subset(dat,cells=which(dat$type.CCC%in%c("BMDM","CD4_naive","CD4_Tcell","Endothelial","Treg")&dat$type_group=="LCBM"))
treg <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
sub.dat$type.CCC[which(colnames(sub.dat)%in%colnames(treg))] <- paste0("Treg_",treg$State)
sub.dat$sample <- sub.dat$orig.ident
sub.dat$sample[grep("T-Bsc1",sub.dat$sample)] <- "T_Bsc1"


#
# BMDM   CD4_naive   CD4_Tcell Endothelial        Treg      Treg_1
# 7433         991         426         231           3         491
# Treg_2
# 429
#=======================================================================================================================

cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}

sample_name <- unique(sub.dat$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/")
    tmp <- subset(sub.dat,cells=which(sub.dat$sample==sample_name[i]))
    cellphoneDB_input(tmp,"type.CCC",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

# run in shell 
for i in `ls /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done

#=======================================================================================================
# check result 
# 2020-12-12
#=======================================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Treg_2.Endothelial")]
	colnames(tmp.value.f) <- c("id_cp_interaction",paste0(sample_name[i],"Treg_2.Endothelial")) 
	cm.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]]$id_cp_interaction%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
# some pair is same id but not same gene 
# ========================================
pair_ann <- tmp.value[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/cm.Treg.pvalue.f.RDS")

#=============================================================================================================================================
# means 
#=============================================================================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.mean.f <- tmp.mean[,c("id_cp_interaction","Treg_2.Endothelial")]
	colnames(tmp.mean.f) <- c("id_cp_interaction",paste0(sample_name[i],"Treg_2.Endothelial")) 
	cm.Treg[[i]] <- tmp.mean.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 

cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]]$id_cp_interaction%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]]$id_cp_interaction%in%co_id),])
}

# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
pair_ann <- tmp.mean[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/cm.Treg.means.f.RDS")

#===================================================================================================================================================
# final calculate result 
# 2020-12-14
#===================================================================================================================================================
pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/cm.Treg.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/cm.Treg.means.f.RDS")
# set some filter 
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]
mean_value <- apply(means.res,1,function(x){mean(x)})
gpair <- names(which(apply(pvalue.res,1,function(x){length(which(x>0))})>1))
mean_value.f <- mean_value[which(names(mean_value)%in%gpair)]



























#==============================================================================
# get cell type refine 
#==============================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/T_cell_featureplot.pdf")
DefaultAssay(dat) <- "RNA"
DimPlot(dat,label=T)
DimPlot(dat,group.by="type_group")
FeaturePlot(dat,features=c("CD4","IL7R","CD3D"),label=T,order=T) # CD4 T
FeaturePlot(dat,features=c("CD8A","CD8B","CD3D"),label=T,order=T) # CD8 T 
FeaturePlot(dat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),label=T,order=T) # Treg 
FeaturePlot(dat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),label=T,order=T) # exhausted
FeaturePlot(dat,features=c("TCF7","SELL","LEF1","CCR7"),label=T,order=T) # naive T
FeaturePlot(dat,features=c("IRF4", "CREM", "NR4A2"),label=T,order=T) # Th17 
dev.off()



#=====================================================================================
# 2020-11-29
#=====================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")

dat$type.rough <- "unsure"
dat$type.rough[which(dat$seurat_clusters%in%c(1,2,5,7,8,9,12))] <- "CD4_Tcell"
dat$type.rough[which(dat$seurat_clusters%in%c(5,8))] <- "Treg"
dat$type.rough[which(dat$seurat_clusters%in%c(1,2,7))] <- "CD4_naive"
dat$type.rough[which(dat$seurat_clusters%in%c(0,3,4,6,10,11))] <- "CD8_Tcell"
dat$type.rough[which(dat$seurat_clusters%in%c(0,4,6,10))] <- "CD8_exhausted/activated"
dat$type.rough[which(dat$seurat_clusters%in%c(13,14))] <- "unsure"

saveRDS(dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")


pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/DimPlot_Tcell.pdf",useDingbats=F)
DimPlot(dat,group.by="type.rough",cols=c("#a0ac48","#8cb811","#ffcc2f","#f1632a","#70b29c","#cecece"))
dev.off()



#=====================================================================================
# sample group different 
# Treg is enrich in LCBM 
#=====================================================================================
table(dat$type.rough,dat$type_group) -> tmp
tmp <- tmp[-6,] # no unsure cells 
apply(tmp,2,function(x){x/sum(x)}) -> tmp1 # 2020-12-17 change 
# apply(tmp1,1,function(x){x/sum(x)}) -> tmp.res 

library(ggplot2)
library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.res <- melt(tmp1)
colnames(tmp.res) <- c("cell_type","type","percentage")
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/type_group_Tcell.pdf",useDingbats=F)
# barplot(tmp.res,ylab = "percentage",xlab = "Cell Type",col = c("#faa918","#7ac70c","#1cb0f6"))
# legend("topright",xpd=T, c("GBM","LC","LCBM"),
# 		pch=c(15,15,15,15,15,15), col = c("#faa918","#7ac70c","#1cb0f6"))
tmp.res$cell_type <- factor(tmp.res$cell_type,levels=rev(c("CD4_naive","CD8_exhausted/activated","CD8_Tcell","CD4_Tcell","Treg")))
ggplot(tmp.res, aes(x = type, y = percentage, fill = cell_type,stratum = cell_type, alluvium = cell_type)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = rev(c("#a0ac48","#ffcc2f","#f1632a","#8cb811","#70b29c")))+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()






#====================================================================================
# use Sangsung daa to verify
#====================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat <- subset(dat,cells=which(dat$Sample_Origin%in%c("mBrain","tLung")&dat$Cell_type=="T lymphocytes"))
table(sub.dat$Cell_subtype,sub.dat$Sample_Origin) -> a
a.f <- a[-which(rowSums(a)==0),]
res <- a.f[,-which(colSums(a.f)==0)]
apply(res,2,function(x){x/sum(x)})
# verify in Sang sung  and could use Barplot to show 


#====================================================================================
# use GSE143423 to verfiy CD4 and ANXA1 T cell 
# FeaturePlot
#====================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE143423/GSE143423_seurat.RDS")
data <- dat[['RNA']]@data
dat$CD4bin <- 0
dat$CD4bin[which(data["CD4",]>0)] <- 1
dat$ANXA1bin <- 0
dat$ANXA1bin[which(data["ANXA1",]>0)] <- 1
#====================================================================================
tmp <- FeaturePlot(dat,features=c("CD4bin","ANXA1bin"),blend=T,combine=F)
g <- tmp[[3]]
g + scale_color_manual(values=c("black","black","black","red"))






#======================================================================================
# Treg development with BMS correlation in TCGA samples 
# lly think when use signature ,just use a treg reactome signature is ok 
# 2020-12-16
#======================================================================================
# use gmt file 
options(stringsAsFactors=F)
dat <- read.delim2("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/cd4Tcell.txt.gmt",sep="\t",header=F)
tmp <- unique(c(as.vector(unlist(dat[8,])),as.vector(unlist(dat[30,])),as.vector(unlist(dat[44,]))))
treg_mod <- tmp[-grep("http|GSE[0-9]+",tmp)]


#================
# treg reactome signature
genes <- c("FOXP3","CBFB","NFATC2","IL2","IFNG","IL2RA","RUNX1","CTLA4","TNFRSF18","CR1")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(genes,'Treg',out='/public/workspace/lily/MOD_file/Treg.mod')

# calculate TCGA LUAD 
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Treg"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)

# load BMS 
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")

#=============================================================================
# plot result 
#=============================================================================
mod$BMS <- luad_mod[,2]
library(ggExtra)
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_BMS_cor.pdf",useDingbats=F)
p<-ggplot(mod,aes(x=BMS,y=Treg_norm)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Treg Signature") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()





#=============================================================================
# use th17 signature to find subtypeI and subtypeII is any difference 
#=============================================================================
dat <- read.delim2("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/cd4Tcell.txt.gmt",sep="\t",header=F)
#tmp <- intersect(as.vector(unlist(dat[18,])),as.vector(unlist(dat[20,])) %>% intersect(as.vector(unlist(dat[22,]))))
tmp <- unique(c(as.vector(unlist(dat[18,])),as.vector(unlist(dat[20,])),as.vector(unlist(dat[22,]))))
th17_mod <- tmp[-grep("http|GSE[0-9]+",tmp)]

# try other signature 
tmp <- (as.vector(unlist(dat[22,])))[-c(1:2)]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(tmp,'th17_nTreg_up',out='/public/workspace/lily/MOD_file/th17_nTreg_up.mod')
#==============================================
# use LCBM Treg cells to ssGSEA
#==============================================
library(Seurat)
tmp_tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tmp_tcell[['RNA']]@data),c("th17_nTreg_up"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)








#=====================================================================================
# 2020-12-17 
# do not subset Treg cell 
# just use all Treg data in LCBM 
#=====================================================================================
library(Seurat)
Tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
sub.dat <- subset(Tcell,cells=which(Tcell$type_group=="LCBM"&Tcell$type.rough=="Treg"))
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
dat$type.tmp <- "others"
dat$type.tmp[which(colnames(dat)%in%colnames(sub.dat))] <- "Treg"

# just find in LCBM samples 
sub.LCBM <- subset(dat,cells=which(dat$type_group=="LCBM"))
sub.LCBM@active.ident <- as.factor(sub.LCBM$type.tmp)
DefaultAssay(sub.LCBM) <- "RNA"
geneset_tmp <- FindMarkers(sub.LCBM,ident.1="Treg",ident.2="others",assay="RNA",logfc.threshold = 0.1)
geneset_tmp <- geneset_tmp[order(geneset_tmp$avg_logFC,decreasing=T),] # maybe you can use logfc >1 p.adj < 0.05 to filter some genes

gene.f <- geneset_tmp[which(geneset_tmp$p_val_adj<0.05&geneset_tmp$avg_logFC>0.5&geneset_tmp$pct.2<0.5),]

# KEGG pathway analysis 
#===================================================================================
genes <- rownames(gene.f)
# KEGG pathway analysis 
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig
geneinfo <- bitr(genes, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")    
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
tmp1 <- ekk@result
tmp1 <- tmp1[which(tmp1$pvalue<0.05),]




#===================================================================================
# run cellphoneDB to find Both TF tumor and Treg 
#===================================================================================

#=====================================================================================
library(Seurat)
Tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
sub.dat <- subset(Tcell,cells=which(Tcell$type_group=="LCBM"&Tcell$type.rough=="Treg"))
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
dat$type.CCC <- as.vector(dat$type)
dat$type.CCC[which(colnames(dat)%in%colnames(sub.dat))] <- "Treg" # treg 
# BMDM and MG 
Mye <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
dat$type.CCC[which(colnames(dat)%in%colnames(Mye))] <- as.vector(Mye$type.refine)

# Both TF 
data <- dat[['RNA']]@data
dat$type.CCC[which(data["MYBL2",]>0&data["CEBPB",]>0&dat$type=="malignant"&dat$type_group=="LCBM")] <- "Both_T"

# get result Seurat and run CellphoneDB 
sub.res <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type.CCC%in%c("BMDM","Both_T","malignant","T_cell","Treg")))

# run CellphoneDB 
#===========================================================================================================================
cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}

sub.res$sample <- sub.res$orig.ident
sub.res$sample[grep("T-Bsc1",sub.res$sample)] <- "T_Bsc1"
sample_name <- unique(sub.res$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/")
    tmp <- subset(sub.res,cells=which(sub.res$sample==sample_name[i]))
    cellphoneDB_input(tmp,"type.CCC",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

# run in shell 
for i in `ls /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done



#====================================================================================================================
# analysis result 
# Treg to Tumor
#====================================================================================================================
# run in other script 















































#=====================================================================================
# immunce cell AI 
# Treg 
# maybe should put in supplmentary
#=====================================================================================
dat <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/ImmuCellAI_abundance_result_TCGA_LUAD.txt",sep="\t",header=T)
colnames(dat)[1] <- "Sample"
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
dat$BMS <- luad_mod[,2]

dat$type <- "unknow"
dat$type[which(dat$BMS>quantile(dat$BMS,0.67))] <- "High"
dat$type[which(dat$BMS<quantile(dat$BMS,0.33))] <- "Low"

apply(dat[,-c(1,28)],2,function(x){wilcox.test(x[which(dat$type=="High")],x[which(dat$type=="Low")])$p.value})



pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TCGA_nTreg_BMS.pdf")
boxplot(nTreg~type,data=dat)
legend("topright",legend=paste0("P=",wilcox.test(nTreg~type,data=dat)$p.value))
dev.off()



#===================================================================================
# immune escape genes calculate mod 
# calculate BMS with immune escape score correlation 
#===================================================================================
dat <- read.csv("/public/workspace/zhumy/ref/immune_evasion_genes.csv")
c(na.omit(unique(as.vector(dat$HGNC.symbol)))) -> gene
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,'immuneEscape',out='/public/workspace/lily/MOD_file/immuneEscape.mod')

# calculate TCGA LUAD 
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("immuneEscape"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)

# result 
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData")
mod$BMS <- luad_mod[,2]

#===================================================================================
library(ggExtra)
library(ggplot2)
library(ggpubr)

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/immuneEscape_BMS_cor.pdf",useDingbats=F)
p<-ggplot(mod,aes(x=BMS,y=immuneEscape_norm)) + 
    stat_density2d(aes(alpha = 1,fill = ..density..),geom = "raster", contour = FALSE) +
    scale_fill_gradient (low = "#FFFFFF", high = "#377EB8") +
    ylab("Immune Escape Signature") + xlab('BMS score') + stat_smooth(method="lm",se=T) + 
    stat_cor(method = "pearson",size=2,label.x = 0.1,label.y =0.1) + 
    geom_point(colour='#377EB8',size=0.5) + theme_bw() +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="none",axis.line = element_line(colour = "black"), axis.text=element_text(colour="black"))  
 
ggMarginal(p, type="density", margins = "both", fill = "#BBDFFB") 
dev.off()






#=====================================================================================
# run monocle for T_cell 
# to check difference
#=====================================================================================
library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(Seurat)
library(monocle)

# CD4T_cell
tmp_dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
dat <- subset(tmp_dat,cells=which(tmp_dat$type.rough%in%c("CD4_Tcell","CD4_naive","Treg")))

# re-inte 
# 2020-12-17 # change 
# ###### Fri Jan 15 20:19:26 CST 2021
#===================================================================================
inte.list <- list() 
samplelist <- unique(dat$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(dat,cells=which(dat$orig.ident==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}

integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=20,k.score = 20,,dims=1:20)
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
inte <- RunUMAP(inte,dims=1:10)
saveRDS(inte,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_trajectory_pre.RDS") #2020-12-17 ###### Fri Jan 15 20:26:17 CST 2021
# #====================================================================================
#######################################################################################################################################################
# ###### Fri Jan 15 20:27:51 CST 2021
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_trajectory_pre.RDS")
data <- as(as.matrix(dat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = dat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#====================================================================================
# 2021-1-13 
# re-inte CD4 T cell can not give a good monocle tree 
# so change use inte all T cell faeature genes 
#====================================================================================
DefaultAssay(dat) <- "integrated"
genes <-  rownames(dat)
#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- genes
cds <- setOrderingFilter(cds, ordering_genes)
# pseudotime 
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)
cds$type.T <- "CD4T"
# saveRDS(cds,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_cds.RDS")

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4_Tcell_trajectory.pdf",useDingbats=F)
plot_cell_trajectory(cds,color_by="type.rough")+scale_colour_manual(values=c("#bff199","#bff199","#70b29c"))+facet_wrap(~type_group, nrow = 3)
cds$ANXA1 <- "non"
cds$ANXA1[which(cds@assayData$exprs["ANXA1",]>0)] <- "yes"
plot_cell_trajectory(cds,color_by="ANXA1")+facet_wrap(~type_group, nrow = 3) # ANXA1 expression 
plot_cell_trajectory(cds,color_by="State")+facet_wrap(~type_group, nrow = 3)
dev.off()


#===================================================================================
# get CD4T cells to analysis difference 
#===================================================================================
cds <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_cds.RDS")
treg <- pData(cds)[which(cds$State%in%c(1,2)&cds$type.rough=="Treg"),c("type.rough","State")]
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
treg.dat <- subset(dat,cells=which(colnames(dat)%in%rownames(treg)&dat$type_group=="LCBM"))
treg.dat$State <- treg$State[which(rownames(treg)%in%colnames(treg.dat))]

# find different genes 
treg.dat@active.ident <- as.factor(treg.dat$State)
DefaultAssay(treg.dat) <- "RNA"
geneset_tmp <- FindMarkers(treg.dat,ident.1=2,ident.2=1,assay="RNA",logfc.threshold = 0.1)
geneset_tmp <- geneset_tmp[order(geneset_tmp$avg_logFC,decreasing=T),] # maybe you can use logfc >1 p.adj < 0.05 to filter some genes



#set filter to make a faeture gene 
# 2021-1-13
####################################################################################
treg.dat@active.ident <- as.factor(treg.dat$State)
DefaultAssay(treg.dat) <- "RNA"
geneset_tmp <- FindMarkers(treg.dat,ident.1=2,ident.2=1,assay="RNA",logfc.threshold = 0.1)
tmp.f <- geneset_tmp[which(geneset_tmp$p_val_adj<0.05),]
tmp.f$pct_FC <- tmp.f$pct.1/tmp.f$pct.2
tmp.f$pct_DF <- tmp.f$pct.1-tmp.f$pct.2
tmp.f[which(tmp.f$avg_logFC>0.5&tmp.f$pct_DF>0.2),] # treg2
tmp2 <- rownames(tmp.f[which(tmp.f$avg_logFC>0.5&tmp.f$pct_DF>0.2),])
write.table(tmp2,file="./tmp2.txt",col.names=F,row.names=F,quote=F)
# then do hallmarks enrich in enrichr and download table

geneset_tmp <- FindMarkers(treg.dat,ident.1=1,ident.2=2,assay="RNA",logfc.threshold = 0.1)
tmp.f <- geneset_tmp[which(geneset_tmp$p_val_adj<0.05),]
tmp.f$pct_FC <- tmp.f$pct.1/tmp.f$pct.2
tmp.f$pct_DF <- tmp.f$pct.1-tmp.f$pct.2

tmp.f[which(tmp.f$avg_logFC>0.5&tmp.f$pct_DF>0.2),] # treg1
tmp1 <- rownames(tmp.f[which(tmp.f$avg_logFC>0.5&tmp.f$pct_DF>0.2),])
write.table(tmp1,file="./tmp1.txt",col.names=F,row.names=F,quote=F)

# just calclulate Hallmrks for all cells 
# not good result 
###### Wed Jan 13 20:07:05 CST 2021
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# dat <- treg.dat[['RNA']]@data
# modlist <- gsub("\\.mod","",dir("/public/workspace/lily/MOD_file/HALLMARK/"))
# mod <- mod.analyze2(as.matrix(dat),modlist,"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
# mod <- data.frame(mod)
# mod$type <- paste0("Subtype",treg.dat$State)
# mod.f <- mod[,c(1:50)]
# ann <- data.frame(row.names=rownames(mod.f),type=paste0("State",treg.dat$State))
# #
# library(pheatmap)
# pheatmap(mod.f)


##################################################################################################################
# plot hallmark result 
###### Wed Jan 13 20:33:06 CST 2021
# barplot for treg2 
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/MSigDB_Hallmark_TregII.txt",sep="\t",header=T)
tmp.dat <- tmp.dat[order(tmp.dat$Odds.Ratio,decreasing=T),]
tmp.dat.f <- tmp.dat[which(tmp.dat$P.value<0.01),]



####### Fri Jan 15 11:19:08 CST 2021
# use ggplot to show 
library(ggplot2)
tmp.dat.f$Term <- factor(tmp.dat.f$Term,levels=as.vector(tmp.dat.f$Term))
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/tregII_circle.pdf")
ggplot(data = tmp.dat.f) + geom_bar(aes(x=as.factor(Term),y=Odds.Ratio), stat="identity", alpha=0.5)+
	theme_classic() + coord_polar() + scale_fill_viridis_c()
dev.off()



pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/tregII.pdf",height=15)
barplot(tmp.dat.f[,c("Odds.Ratio")],names=tmp.dat.f$Term,las=2)
dev.off()

# barplot for treg1
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/MSigDB_Hallmark_TregI.txt",sep="\t",header=T)
tmp.dat <- tmp.dat[order(tmp.dat$Odds.Ratio,decreasing=T),]
tmp.dat.f <- tmp.dat[which(tmp.dat$P.value<0.01),]




library(ggplot2)
tmp.dat.f$Term <- factor(tmp.dat.f$Term,levels=as.vector(tmp.dat.f$Term))
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/tregI_circle.pdf")
ggplot(data = tmp.dat.f) + geom_bar(aes(x=as.factor(Term),y=Odds.Ratio), stat="identity", alpha=0.5)+
	theme_classic() + coord_polar() + scale_fill_viridis_c()
dev.off()




pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/tregI.pdf",height=15)
barplot(tmp.dat.f[,c("Odds.Ratio")],names=tmp.dat.f$Term,las=2)
dev.off()




###################################################################################################
# T cell trans ednothelial
###### Thu Jan 14 10:07:37 CST 2021
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
gene <- c("ACTB","ACTG1","ACTN1","ACTN2","ACTN3","ACTN4",
	"AFDN","ARHGAP35","ARHGAP5","BCAR1","CD99","CDC42","CDH5",
	"CLDN1","CLDN10","CLDN11","CLDN14","CLDN15","CLDN16",
	"CLDN17","CLDN18","CLDN19","CLDN2","CLDN20","CLDN22",
	"CLDN23","CLDN3","CLDN4","CLDN5","CLDN6","CLDN7","CLDN8",
	"CLDN9","CTNNA1","CTNNA2","CTNNA3","CTNNB1","CTNND1",
	"CXCL12","CXCR4","CYBA","CYBB","ESAM","EZR","F11R","GNAI1",
	"GNAI2","GNAI3","ICAM1","ITGA4","ITGAL","ITGAM","ITGB1",
	"ITGB2","ITK","JAM2","JAM3","MAPK11","MAPK12","MAPK13","MAPK14",
	"MMP2","MMP9","MSN","MYL10","MYL12A","MYL12B","MYL2","MYL5",
	"MYL7","MYL9","MYLPF","NCF1","NCF2","NCF4","NOX1","NOX3","OCLN",
	"PECAM1","PIK3CA","PIK3CB","PIK3CD","PIK3CG","PIK3R1","PIK3R2",
	"PIK3R3","PIK3R5","PLCG1","PLCG2","PRKCA","PRKCB","PRKCG","PTK2",
	"PTK2B","PTPN11","PXN","RAC1","RAC2","RAP1A","RAP1B","RAPGEF3",
	"RAPGEF4","RASSF5","RHOA","RHOH","ROCK1","ROCK2","SIPA1","THY1",
	"TXK","VASP","VAV1","VAV2","VAV3","VCAM1","VCL"
)
mod.generate(gene,'TcellTransEndo',out='/public/workspace/lily/MOD_file/TcellTransEndo.mod')

##################################################################################################
cds <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_cds.RDS")
treg <- pData(cds)[which(cds$State%in%c(1,2)&cds$type.rough=="Treg"),c("type.rough","State")]
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
treg.dat <- subset(dat,cells=which(colnames(dat)%in%rownames(treg)&dat$type_group=="LCBM"))
treg.dat$State <- treg$State[which(rownames(treg)%in%colnames(treg.dat))]
saveRDS(treg.dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
#=================================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
mod <- mod.analyze2(as.matrix(dat[['RNA']]@data),c("TcellTransEndo"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- data.frame(mod)
mod$type <- paste0("State",dat$State)

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TtransEndo.pdf",useDingbats=F)
boxplot(TcellTransEndo_norm~type,data=mod,FUN=median,outliers=F)
dev.off()










############################################################################################################################################
# ###### Fri Jan 15 19:09:06 CST 2021
# GSE131907 verify
# check Treg
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
tmp.sub <- subset(dat,cells=which(dat$Cell_subtype=="Treg"))
sub.dat <- subset(tmp.sub,cells=which(tmp.sub$Sample_Origin=="mBrain"|tmp.sub$Sample_Origin=="tLung"))

# re-cluster 
recluster <- function(tmp_dat){
# seurat object
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
	# PCA
	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:20)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=1)
	# T-SNE
	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:20,check_duplicates = FALSE)
	return(tmp_dat)
}
sub.dat <- recluster(sub.dat)
































#===================================================================================
# save result 
# plot violin plot 
#===================================================================================

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4_Tcell_ANXAN1.pdf",useDingbats=F)
VlnPlot(treg.dat,features="ANXA1")
dev.off()

saveRDS(treg.dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")



# try to find marker gene in CD4_ANXA1 T cell 
#===================================================================================
library(Seurat)
tmp_tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
dat$type.tmp <- "others"
dat$type.tmp[which(colnames(dat)%in%colnames(tmp_tcell)[which(tmp_tcell$State==2)])] <- "ANXA1_CD4T"

# just find in LCBM samples 
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))

sub.dat@active.ident <- as.factor(sub.dat$type.tmp)
DefaultAssay(sub.dat) <- "RNA"
geneset_tmp <- FindMarkers(sub.dat,ident.1="ANXA1_CD4T",ident.2="others",assay="RNA",logfc.threshold = 0.1)
geneset_tmp <- geneset_tmp[order(geneset_tmp$avg_logFC,decreasing=T),] # maybe you can use logfc >1 p.adj < 0.05 to filter some genes

gene.f <- geneset_tmp[which(geneset_tmp$p_val_adj<0.05&geneset_tmp$avg_logFC>0.5&geneset_tmp$pct.1>0.5&geneset_tmp$pct.2<0.5),]






# 2020-12-17
# use another type treg to find marker gene 
#===============================================================================================
library(Seurat)
tmp_tcell <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
dat$type.tmp <- "others"
dat$type.tmp[which(colnames(dat)%in%colnames(tmp_tcell)[which(tmp_tcell$State==1)])] <- "CD4T_1"

# just find in LCBM samples 
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))

sub.dat@active.ident <- as.factor(sub.dat$type.tmp)
DefaultAssay(sub.dat) <- "RNA"
geneset_tmp <- FindMarkers(sub.dat,ident.1="CD4T_1",ident.2="others",assay="RNA",logfc.threshold = 0.1)
geneset_tmp <- geneset_tmp[order(geneset_tmp$avg_logFC,decreasing=T),] # maybe you can use logfc >1 p.adj < 0.05 to filter some genes

gene.f <- geneset_tmp[which(geneset_tmp$p_val_adj<0.05&geneset_tmp$avg_logFC>0.5&geneset_tmp$pct.2<0.5),]









# make a siganture maybe 
#===================================================================================
genes <- rownames(gene.f)
# KEGG pathway analysis 
library(clusterProfiler)
library(org.Hs.eg.db)
require(DOSE)
require(doseplot)
library(igraph)
library(ggplot2)
# sig
geneinfo <- bitr(genes, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")    
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
tmp1 <- ekk@result
tmp1 <- tmp1[which(tmp1$pvalue<0.05),]

#=================================================================================
# use these genes to generate a mod 
# then calculate with angiogenes correlation 
#=================================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(genes,'ANXATcell',out='/public/workspace/lily/MOD_file/ANXA1Tcell.mod')

# calculate TCGA LUAD 
load("/public/workspace/lily/TCGA_Data/TCGA_LUAD/TCGA_LUAD_data.RData")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("ANXA1Tcell","blood_signature"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)


#================================================================================
# do GSEA use  in local 
#================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.RDS")
# make a gct file 
tmp <- as.matrix(dat[["RNA"]]@data)
tmp.res <- cbind(data.frame(NAME=rownames(tmp),Description=rep("NA",nrow(tmp))),tmp)
write.table(tmp.res,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.gct",sep="\t",quote=F,row.names=F)

# sed -i '1 i 25347\t920' /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.gct
# sed -i '1 i #1.2' /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.gct


# make a gls file 
cat(c(920,2,1),file="~/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.cls",sep="\t")
cat("\n",file="~/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.cls",append=T)
cat(c("#State1","State2"),file="~/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.cls",sep="\t",append=T)
cat("\n",file="~/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.cls",append=T)
cat(as.numeric(dat$State),file="~/Lung2Brain/TME/Final_11_28/T_cell/Treg_LCBM.cls",sep="\t",append=T)








#=========================================================================================================
# cellphoneDB result in BMDM with T cells 
# just in LCBM samples
# 2020-12-11
#========================================================================================================

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/")
BMDM.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","BMDM.Treg")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".BMDM.Treg") 
	BMDM.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(BMDM.Treg)){
	tmp.id <- c(tmp.id,as.vector(BMDM.Treg[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



BMDM.Treg.res <- cbind(BMDM.Treg[[1]][which(BMDM.Treg[[1]]$id_cp_interaction%in%co_id),],
							BMDM.Treg[[2]][which(BMDM.Treg[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(BMDM.Treg)){
	BMDM.Treg.res <- cbind(BMDM.Treg.res,BMDM.Treg[[i]][which(BMDM.Treg[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
BMDM.Treg.res <- BMDM.Treg.res[,c(1,grep("BMDM",colnames(BMDM.Treg.res)))]
# add pair gene information 
Communcation.res <- merge(BMDM.Treg.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/BMDM.Treg.pvalue.f.RDS")





#===========================================================================================================================
# GSE 36176 to verify signature and anti angiogenesis therapy 
#===========================================================================================================================
dat <- read.table("~/metastasis/data/verify/GSE36176/GSE36176_series_matrix.txt",sep="\t",header=T,comment.char="!")

# trans id 
gpl570 <- read.table("/public/workspace/lily/metastasis/data/verify/GPL570.txt",sep="\t",header=T)
tmp.res <- merge(dat,gpl570,by.x="ID_REF",by.y="ID")
tmp.res.f <- tmp.res[-grep("///",tmp.res$Gene.Symbol),]
tmp.res.f$ID_REF <- NULL
res <- aggregate(.~Gene.Symbol,data=tmp.res.f,FUN=median)
res.final <- res[-c(1:25),]
rownames(res.final) <- res.final$Gene.Symbol
res.final$Gene.Symbol <- NULL
saveRDS(res.final,file="~/metastasis/data/verify/GSE36176/GSE36176_expr.RDS")

#============================================================================
# calculate ANXA1 signature 
# 2020-12-15
#============================================================================
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(res.final),c("ANXA1Tcell","blood_signature"),"/public/workspace/lily/MOD_file/",permN=1000)
mod <- data.frame(mod)
mod$sample <- rownames(mod)
ann <- read.table("~/metastasis/data/verify/GSE36176/GSE36176_anno.txt",sep="\t") # sample info 
colnames(ann) <- c("sample","type")
merge.res <- merge(ann,mod,by="sample")




















#===================================================================================================================================================
# CD8T_cell
#===================================================================================================================================================
tmp_dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
dat <- subset(tmp_dat,cells=which(tmp_dat$type.rough%in%c("CD8_Tcell","CD8_exhausted")))


#====================================================================================
data <- as(as.matrix(dat@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = dat@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#====================================================================================
# get order genes 
#
#====================================================================================
DefaultAssay(dat) <- "integrated"
genes <-  rownames(dat)
#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

ordering_genes <- genes
cds <- setOrderingFilter(cds, ordering_genes)
# pseudotime 
cds <- reduceDimension(cds,method = 'DDRTree')
cds <- orderCells(cds)
cds$type.T <- "CD8T"
saveRDS(cds,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD8T_cds.RDS")






#=======================================================================================
# plot result 
# CD4T cell 
#=======================================================================================
cds <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_cds.RDS")
plot_cell_trajectory(cds,color_by="type_group")
cds$type.refine <- "unsure"
cds$type.refine[which(cds$seurat_clusters%in%c(5,8))] <- "Treg"
cds$type.refine[which(cds$seurat_clusters%in%c(1,2,7))] <- "naive_CD4"
plot_cell_trajectory(cds,color_by="type.refine")

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CD4T_trajectoy.pdf")
plot_cell_trajectory(cds,color_by="type_group")
plot_cell_trajectory(cds,color_by="type.refine")+scale_colour_manual(values=c("#76b852","#ee3322","#eff3f6"))+facet_wrap(~type_group, nrow = 3)
dev.off()


#===============================
# treg ?
#===============================
treg <- pData(cds)[which(cds$State%in%c(1,2)&cds$type.refine=="Treg"),c("type.refine","State")]

dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
treg.dat <- subset(dat,cells=which(colnames(dat)%in%rownames(treg)&dat$type_group=="LCBM"))
treg.dat$State <- treg$State[which(rownames(treg)%in%colnames(treg.dat))]


treg.dat@active.ident <- as.factor(treg.dat$State)
geneset_tmp <- FindMarkers(treg.dat,ident.1=2,ident.2=1,assay="RNA",logfc.threshold = 0.1)
geneset_tmp <- geneset_tmp[order(geneset_tmp$avg_logFC,decreasing=T),]

#===============================
# # a try 
#===============================
tmp.data <- treg.dat[['RNA']]@data
treg.dat$marker1 <- "no"
treg.dat$marker1[which(tmp.data["CD4",]>0&tmp.data["ANXA1",]>0)] <- "yes"
table(treg.dat$marker1)

# pdf("./tmp.pdf")
# DimPlot(treg.dat,group.by="State")
# DefaultAssay(treg.dat) <- "RNA"
# FeaturePlot(treg.dat,features=c("IRF4", "CREM", "NR4A2"))
# dev.off()

# subset.dat <- subset(treg.dat,cells=which(treg.dat$State==2))
tmp.data <- subset.dat[['RNA']]@data
subset.dat$marker1 <- "no"
subset.dat$marker1[which(tmp.data["CTLA4",]>0&tmp.data["KRT18",]>0)] <- "yes"
table(subset.dat$marker1)









#======================================================================================================
# BMDM 2 Treg 
# 2020-12-31
#======================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
BMDM.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","BMDM.Treg")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".BMDM.Treg") 
	BMDM.Treg[[i]] <- tmp.value.f
}
#######################################################################################################
tmp.id <- c()
for(i in 1:length(BMDM.Treg)){
	tmp.id <- c(tmp.id,as.vector(BMDM.Treg[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 
#######################################################################################################
BMDM.Treg.res <- cbind(BMDM.Treg[[1]][which(BMDM.Treg[[1]]$id_cp_interaction%in%co_id),],
							BMDM.Treg[[2]][which(BMDM.Treg[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(BMDM.Treg)){
	BMDM.Treg.res <- cbind(BMDM.Treg.res,BMDM.Treg[[i]][which(BMDM.Treg[[i]]$id_cp_interaction%in%co_id),])
}

# do some modify
BMDM.Treg.res <- BMDM.Treg.res[,c(1,grep("Treg",colnames(BMDM.Treg.res)))]
# add pair gene information 
Communcation.res <- merge(BMDM.Treg.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/BMDM_Treg_pvalue.f.RDS")

#################################################################
# means result 
# 2020-12-31
#================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
BMDM.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","BMDM.Treg")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".BMDM.Treg") 
	BMDM.Treg[[i]] <- tmp.value.f
}
#######################################################################################################
tmp.id <- c()
for(i in 1:length(BMDM.Treg)){
	tmp.id <- c(tmp.id,as.vector(BMDM.Treg[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 
#######################################################################################################
BMDM.Treg.res <- cbind(BMDM.Treg[[1]][which(BMDM.Treg[[1]]$id_cp_interaction%in%co_id),],
							BMDM.Treg[[2]][which(BMDM.Treg[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(BMDM.Treg)){
	BMDM.Treg.res <- cbind(BMDM.Treg.res,BMDM.Treg[[i]][which(BMDM.Treg[[i]]$id_cp_interaction%in%co_id),])
}

# do some modify
BMDM.Treg.res <- BMDM.Treg.res[,c(1,grep("Treg",colnames(BMDM.Treg.res)))]
# add pair gene information 
Communcation.res <- merge(BMDM.Treg.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/BMDM_Treg_mean.f.RDS")

######################################################################################################################################
# final check result 
# 2020-12-31
######################################################################################################################################
pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/BMDM_Treg_pvalue.f.RDS")
means <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/BMDM_Treg_mean.f.RDS")

means.f <- means[which(rownames(means)%in%rownames(pvalue)),]
which.max(apply(means.f,1,function(x){mean(x)}))
apply(means.f,1,function(x){mean(x)})[order(apply(means.f,1,function(x){mean(x)}))]


#=====================================================
# plot result 
# 2021-1-3
######################################################
library(ggplot2)
library(complexpheatmap)
#=====================================================
plot_Circ_heatmap <- function(mat, cluster) {
    # ref: https://zhuanlan.zhihu.com/p/136138642
    # @mat: row(sample or groups) X col(pair info or genes)
    # @cluster: set cluster number in hclust function
    library(dendextend)
    library("circlize")
    library(RColorBrewer)

    mat=scale(mat, center = TRUE, scale = TRUE)
    dend <-as.dendrogram(hclust(dist(t(mat))))
    n=1
    dend <-dend %>% set("branches_k_color", k = n) 
    # par(mar=c(7.5,3,1,0))
    # plot(dend)   # 聚类后的样本信息
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
        canvas.xlim =c(-1,1),
        canvas.ylim = c(-1,1),
        cell.padding = c(0,0,0,0), 
        gap.degree =90
    )
    factors = "a"
    circos.initialize(factors, xlim = c(0, ncol(mat2)+1))
    circos.track(
        ylim = c(0, nr),bg.border = NA,track.height = 0.01*nr,
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
                    cex = 0.1,
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
    # 添加树
    max_height <-max(attr(dend, "height"))
    circos.track(ylim = c(0, max_height),bg.border = NA,track.height = 0.1, 
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


















################################################################################################################
# 2021-1-9
# Treg analysis
#===============================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")

# Run library(dorothea) find which TF activated 
################################################################################################################
# R 4.0.2
library(dorothea)
library(Seurat)
library(bcellViper)
library(dplyr)
library(viper)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
dat <- as.matrix(tmp.dat[["RNA"]]@data)

data(dorothea_hs,package="dorothea")
#subset DoRothEA to the confidence levels A and B to include only the high quality regulons
regulons = dorothea_hs[which(dorothea_hs$confidence%in%c("A","B","C","D")),]
tf_activities <- run_viper(dat, regulons, 
                           options =  list(method = "scale", minsize = 4, 
                                           eset.filter = FALSE, cores = 1, 
                                           verbose = FALSE))




####################################################
# find expression genes to plot heatmap
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
DefaultAssay(dat) <- "RNA"
dat@active.ident <- factor(dat$type.rough)
tmp.res <- FindAllMarkers(dat,assay="RNA")
head(tmp.res[which(tmp.res$cluster=="Treg"&tmp.res$avg_logFC>0),])





###########################################################
# 2021-1-11
# 
# GSE131907
###########################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/GSE131907_all_cell_v12_3.RDS")
sub.dat <- subset(dat,cells=which(dat$Sample_Origin=="mBrain"&dat$Cell_type.refined=="T/NK cells"))

# calculate th17/ angiogenesis features 
dat <- as.matrix(sub.dat[['RNA']]@data)
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("Treg","blood_signature","th17_iTreg_up"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- data.frame(mod)











#################################################################################################
# 2021-1-12
# GSE161116 test 
#================================================================================================
tmp.dat <- read.table("/public/workspace/lily/metastasis/data/verify/GSE161116/GSE161116_series_matrix.txt",sep="\t",comment.char="!",header=T)
ann <- t(read.table("/public/workspace/lily/metastasis/data/verify/GSE161116/sampleinfo.txt"))
tmp.ann <- gsub(" ",".",as.vector(gsub("Patient ","P",ann[2:29,1])))
rownames(tmp.dat) <- tmp.dat$ID_REF
tmp.dat$ID_REF <- NULL
colnames(tmp.dat) <- tmp.ann


# 1. check BMS gene expression 
gene <- readRDS("/public/workspace/lily/Lung2Brain/inte7/BMS_gene.RDS")
tmp.dat[which(rownames(tmp.dat)%in%gene),] -> tmp.1
apply(tmp.1,1,function(x){
	c(	median(x[grep("brain",colnames(tmp.dat))]),
		median(x[grep("lung",colnames(tmp.dat))]),
		wilcox.test(x[grep("brain",colnames(tmp.dat))],x[grep("lung",colnames(tmp.dat))])$p.value
	)
})


# 2. Treg gene expression
tmp.dat[c("IL2RA","FOXP3","IL17A"),] -> tmp.1
apply(tmp.1,1,function(x){
	c(	median(x[grep("brain",colnames(tmp.dat))]),
		median(x[grep("lung",colnames(tmp.dat))]),
		wilcox.test(as.numeric(x[grep("brain",colnames(tmp.dat))]),as.numeric(x[grep("lung",colnames(tmp.dat))]))$p.value
	)
})


# 3. DEG analysis 
tmp.res <-data.frame(t(apply(tmp.dat,1,function(x){
	c(	median(x[grep("brain",colnames(tmp.dat))]),
		median(x[grep("lung",colnames(tmp.dat))]),
		wilcox.test(as.numeric(x[grep("brain",colnames(tmp.dat))]),as.numeric(x[grep("lung",colnames(tmp.dat))]))$p.value
	)
})))
tmp.res$logFC <- log2(tmp.res$brain/tmp.res$lung)
colnames(tmp.res) <- c("brain","lung","p.value")




























































############################################################################################################################################
###### Fri Jan 15 21:33:32 CST 2021
# 2021-1-15 use SingleR to run Cell type identifiy
#===========================================================================================================================================
# library(SingleR)
# library(Seurat)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")

# mid.se <- MonacoImmuneData()
# Tcells_sr <- SingleR(test = as.matrix(dat@assays$RNA@data), ref = mid.se,cluster=seurat_clusters,labels = mid.se$label.fine)
# dat$monacoT <- Tcells_sr$pruned.labels
# # another data to use 
# dice.se <- DatabaseImmuneCellExpressionData()
# Tcells_sr2 <- SingleR(test = as.matrix(dat@assays$RNA@data), ref = dice.se, 
#     labels = dice.se$label.fine)
# dat$dice <- Tcells_sr2$pruned.labels

# #==========================================================================================================================================
# # SingleR
# # re-cluster make a bigger resolution 
# library(Seurat)
# library(SingleR)
# dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
# recluster <- function(tmp_dat){
# # seurat object
# 	tmp_dat <- FindVariableFeatures(object = tmp_dat)
# 	# scaling
# 	all.genes <- rownames(x = tmp_dat)
# 	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
# 	# PCA
# 	tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
# 	# clustering
# 	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:20)
# 	# select proper resolution
# 	tmp_dat <- FindClusters(object = tmp_dat,resolution=2)
# 	# T-SNE
# 	tmp_dat <- RunTSNE(object = tmp_dat,dims=1:20,check_duplicates = FALSE)
# 	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
# 	return(tmp_dat)
# }
# dat <- recluster(dat)
# mid.se <- MonacoImmuneData()
# Tcells_sr <- SingleR(test = as.matrix(dat@assays$RNA@data), ref = mid.se,cluster=seurat_clusters,labels = mid.se$label.fine)
# dat$monacoT <- Tcells_sr$pruned.labels
# # check percentage 
# tmp <- table(dat$seurat_clusters,dat$monacoT)
# apply(tmp,1,function(x){x/sum(x)})








########################################################################################################################################
# Single R needs a best refernce 
# so just find marker and use marker to identify cluster
########################################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
Deg <- FindAllMarkers(dat,assay="RNA",only.pos=T) # Find DEG 
Deg.f <- Deg[which(Deg$p_val_adj<0.05),]
Deg.f <- Deg.f[order(Deg.f$avg_logFC,decreasing=T),]

# ###### Mon Jan 18 11:10:13 CST 2021
# Finally use cell type 
#=======================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
dat$celltype.refine <- "Unclassify"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==0)] <- "activated/exhausted CD8"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==1)] <- "naive CD4"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==2)] <- "exhausted CD4"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==3)] <- "CD4/CD8 mix"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==4)] <- "cytotoxic CD8"
dat$celltype.refine[which(dat$integrated_snn_res.0.8%in%c(5,8))] <- "Treg"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==6)] <- "proliferation CD8"
dat$celltype.refine[which(dat$integrated_snn_res.0.8==7)] <- "follicular Th"
dat$celltype.refine[which(dat$integrated_snn_res.0.8%in%c(10))] <- "cytotoxic/proliferation CD8"
dat$celltype.refine[which(dat$integrated_snn_res.0.8%in%c(11,12))] <- "NK"

###########################################################################
#   activated/exhausted CD8                 CD4/CD8 mix
#                        1871                        1360
#               cytotoxic CD8 cytotoxic/proliferation CD8
#                        1261                         450
#               exhausted CD4               follicular Th
#                        1361                         603
#                   naive CD4                          NK
#                        1423                         756
#           proliferation CD8                        Treg
#                         609                        1610
#                  Unclassify
#                         616

cols <- c("#85b7e2","#ef9020","#91be3e","#e990ab","#ffd616","#5c92fa","#589636","#1fadc5","#9013fe","#ec4534","#969491")
DimPlot(dat,group.by="celltype.refine",cols=cols)


#################################################################################################################################
#===========================================================================
# ###### Mon Jan 18 14:46:41 CST 2021
#===========================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
# 1. DimPlot result 
cols <- c("#85b7e2","#ef9020","#91be3e","#e990ab","#ffd616","#5c92fa","#589636","#1fadc5","#9013fe","#ec4534","#969491")
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/DimPlot_Tcell.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype.refine",cols=cols)
dev.off()

# 2. cell type percentage change 
table(dat$celltype.refine,dat$type_group) -> tmp
tmp <- tmp[-11,] # no unclassify cells 
apply(tmp,2,function(x){x/sum(x)}) -> tmp1 # 2020-12-17 change 
# apply(tmp1,1,function(x){x/sum(x)}) -> tmp.res 

library(ggplot2)
library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.res <- melt(tmp1)
colnames(tmp.res) <- c("cell_type","type","percentage")
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/type_group_Tcell.pdf",useDingbats=F)
# barplot(tmp.res,ylab = "percentage",xlab = "Cell Type",col = c("#faa918","#7ac70c","#1cb0f6"))
# legend("topright",xpd=T, c("GBM","LC","LCBM"),
# 		pch=c(15,15,15,15,15,15), col = c("#faa918","#7ac70c","#1cb0f6"))
# tmp.res$cell_type <- factor(tmp.res$cell_type,levels=rev(c("CD4_naive","CD8_exhausted/activated","CD8_Tcell","CD4_Tcell","Treg")))
ggplot(tmp.res, aes(x = type, y = percentage, fill = cell_type,stratum = cell_type, alluvium = cell_type)) +
geom_stratum() +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()





#############################################################################################################
# CD4 T cell subset run monocle 
#============================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
sub.dat <- subset(dat,cells=which(dat$celltype.refine%in%c("exhausted CD4","follicular Th","naive CD4","Treg")))

# run monocle 
# 1. re-inte
inte.list <- list() 
samplelist <- unique(sub.dat$orig.ident)
for(i in 1:length(samplelist)){
	tmp <- subset(sub.dat,cells=which(sub.dat$orig.ident==samplelist[i]))
	DefaultAssay(tmp) <- "RNA"
	inte.list[i] <- tmp
}

# Find variable genes
integration.anchors <- FindIntegrationAnchors(object.list =inte.list,dims=1:10,k.filter=10,k.score=10)
inte <- IntegrateData(anchorset = integration.anchors,dims=1:10,k.weight=10)
# inte <- FindVariableFeatures(inte)
order_gene <- rownames(inte)
#############################################################################################################
# run monocle 
library(monocle)
DefaultAssay(sub.dat) <- "RNA"
cds.dat <- Seurat::as.CellDataSet(sub.dat)
cds.dat <- estimateSizeFactors(cds.dat)
cds.dat <- estimateDispersions(cds.dat)
cds.dat <- detectGenes(cds.dat, min_expr = 0.1)
cds.dat <- setOrderingFilter(cds.dat, ordering_genes=order_gene)
# pseudotime 
cds <- reduceDimension(cds.dat,method = 'DDRTree')
cds <- orderCells(cds)


plot_cell_trajectory(cds,color_by="celltype.refine")+scale_colour_manual(values=c("#ffd616","#5c92fa","#589636","#ec4534"))+facet_wrap(~type_group, nrow = 3)



# try to run Hallmarks enrichment analysis 
sub.dat$State <- cds$State
tmp <- subset(sub.dat,cells=which(sub.dat$celltype.refine=="Treg"&sub.dat$type_group=="LCBM"))
gene <- FindMarkers(tmp,ident.1=1,ident.2=3,assay="RNA")
gene.f <- gene[which(gene$p_val_adj<0.05&gene$pct.2<0.5),]
gene.f[which(gene.f$avg_logFC>1),]
write.table(rownames(gene.f[which(gene.f$avg_logFC>1),]),file="./tmp_gene.txt",quote=F,row.names=F,col.names=F)

# calculate T transe endothelials
# 2021-1-18
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(tmp[['RNA']]@data),c("TcellTransEndo"),"/public/workspace/lily/MOD_file/",permN=0)








































































#=====================================================================================================================================================
# ###### Thu Jan 21 19:51:37 CST 2021
# a new way to anallysis T cell 
######################################################################################################################################################
# 1. calculate cell cell communication
# Treg and BMDM/MG
# you can find this code in Treg_CCC code file by find "2021-1-21"

# check Final result 
# first is BMDM
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]



# Venne Plot of result 
library(venn)
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Myeloid_Treg_Venn.pdf")
venn(list(rownames(MDM.mean),rownames(MG.mean)),zcolor = 'style')
dev.off()

# intersection gene pair 
# remove ab complex
inte <- intersect(rownames(MDM.mean),rownames(MG.mean))
tmp.gene <- as.vector(sapply(strsplit(sapply(strsplit(inte,"\\."),function(x){x[2]}),"_"),function(x){c(x[1],x[2])}))
write.table(unique(tmp.gene),file="./tmpgene.txt",row.names=F,col.names=F,quote=F)



# MDM specific gene 
tmp <- rownames(MDM.mean[-which(rownames(MDM.mean)%in%intersect(rownames(MDM.mean),rownames(MG.mean))),])
tmp.gene <- as.vector(sapply(strsplit(sapply(strsplit(tmp,"\\."),function(x){x[2]}),"_"),function(x){c(x[1],x[2])}))
write.table(unique(tmp.gene),file="./tmpgene.txt",row.names=F,col.names=F,quote=F)



#=========================================================================================================================
# plot circle heatmap
# ###### Sat Jan 23 17:28:54 CST 2021
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]




#=======================================================================================================================
# plot sangji 
# ###### Mon Jan 25 10:23:17 CST 2021
#=======================================================================================================================
library(networkD3)
# first is BMDM
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]

# MDM specific interation pairs 
MDM.sp <- MDM.mean[-which(rownames(MDM.mean)%in%intersect(rownames(MDM.mean),rownames(MG.mean))),]
# change format to plot 
tmp.dat <- as.data.frame(apply(MDM.sp,1,function(x){mean(x)}))
colnames(tmp.dat) <- "value"
tmp.dat$source <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[1]})
tmp.dat$target <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[2]})
rownames(tmp.dat) <- NULL
sankeynode <- data.frame(row.names=unique(c(tmp.dat$source,tmp.dat$target)),name=unique(c(tmp.dat$source,tmp.dat$target)))
rownames(sankeynode) <- NULL
sankeynode$ID <- 0:(nrow(sankeynode)-1)

# merge data to get ID 
tmp.1 <- merge(tmp.dat,sankeynode,by.x="source",by.y="name")
colnames(tmp.1)[4] <- "from"
tmp.2 <- merge(tmp.2,sankeynode,by.x="target",by.y="name")
colnames(tmp.2)[5] <- "to"

# final get data 
sankeydata <- tmp.2[,c("from","to","value")]
rownames(sankeynode) <- sankeynode$ID

p <- sankeyNetwork(Links = sankeydata, Nodes = sankeynode, Source = "from",
              Target = "to", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30
)






# Sankey plot is not very good , try to use other plot 
# ###### Mon Jan 25 15:01:30 CST 2021
#======================================================================================================================================================
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]

# MDM specific interation pairs 
MDM.sp <- MDM.mean[-which(rownames(MDM.mean)%in%intersect(rownames(MDM.mean),rownames(MG.mean))),]


# connect 
tmp.dat <- as.data.frame(apply(MDM.sp,1,function(x){mean(x)}))
colnames(tmp.dat) <- "value"
tmp.dat$source <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[1]})
tmp.dat$target <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[2]})
rownames(tmp.dat) <- NULL

# use circlize to plot 
tmp.mat <- matrix(0,nrow=length(unique(tmp.dat$source)),ncol=length(unique(tmp.dat$target)))
rownames(tmp.mat) <- unique(tmp.dat$source)
colnames(tmp.mat) <- unique(tmp.dat$target)

for(i in 1:nrow(tmp.dat)){
	s <- tmp.dat[i,"source"]
	t <- tmp.dat[i,"target"]
	tmp.mat[s,t] <- 1
}

library(circlize)
library(dplyr)
library(readr)
library(stringr)

rece <- rep("grey",ncol(mat))
names(rece) <- colnames(mat)
grid.col = c(LGALS9 = "red",rece) # set color

tmp.1 = rep("ligand",nrow(mat))  # set group
names(tmp.1) <- rownames(mat)
tmp.2 = rep("receproter",ncol(mat))  # set group
names(tmp.2) <- colnames(mat)
group <- c(tmp.1,tmp.2)

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/MDM_Treg_sp_circle.pdf")
chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"),annotationTrackHeight = c(0.03, 0.01),
	directional = 1,direction.type= "arrows",transparency = 0.5,big.gap = 30, small.gap = 1	,group=group
)
circos.clear()
dev.off()









#==================================================================================================================
# MG specific 
#==================================================================================================================
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_Treg.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_Treg.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]

# MDM specific interation pairs 
MG.sp <- MG.mean[-which(rownames(MG.mean)%in%intersect(rownames(MDM.mean),rownames(MG.mean))),]

# change format to a dataframe
tmp.dat <- as.data.frame(apply(MG.sp,1,function(x){mean(x)}))
colnames(tmp.dat) <- "value"
tmp.dat$source <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[1]})
tmp.dat$target <- sapply(strsplit(sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]}),"_"),function(x){x[2]})
rownames(tmp.dat) <- NULL

# use circlize to plot 
tmp.mat <- matrix(0,nrow=length(unique(tmp.dat$source)),ncol=length(unique(tmp.dat$target)))
rownames(tmp.mat) <- unique(tmp.dat$source)
colnames(tmp.mat) <- unique(tmp.dat$target)
# make a matrix 
for(i in 1:nrow(tmp.dat)){
	s <- tmp.dat[i,"source"]
	t <- tmp.dat[i,"target"]
	tmp.mat[s,t] <- 1
}

###############################################################################################################
library(circlize)
library(dplyr)
library(readr)
library(stringr)

mat <- tmp.mat
rece <- rep("grey",ncol(mat))
names(rece) <- colnames(mat)
grid.col = c(LGALS9 = "red",rece) # set color

tmp.1 = rep("ligand",nrow(mat))  # set group
names(tmp.1) <- rownames(mat)
tmp.2 = rep("receproter",ncol(mat))  # set group
names(tmp.2) <- colnames(mat)
group <- c(tmp.1,tmp.2)


pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/MG_Treg_sp_circle.pdf")
chordDiagram(mat, grid.col = grid.col, annotationTrack = c("name", "grid"),annotationTrackHeight = c(0.03, 0.01),
	directional = 1,direction.type= "arrows",transparency = 0.5,big.gap = 30, small.gap = 1	,group=group
)
circos.clear()
dev.off()


#=====================================================================================================================
# add a naive T cell analysis 
# 2021-1-26
#=====================================================================================================================
MDM.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_naiveCD4.mean.f.RDS")
MDM.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/BMDM_naiveCD4.pvalue.f.RDS")
MDM.mean <- MDM.tmp.mean[which(rownames(MDM.tmp.mean)%in%rownames(MDM.pvalue)),]
# other is MG
MG.tmp.mean <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_naiveCD4.mean.f.RDS")
MG.pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_naiveCD4.pvalue.f.RDS")
MG.mean <- MG.tmp.mean[which(rownames(MG.tmp.mean)%in%rownames(MG.pvalue)),]

# MDM specific interation pairs 
MDM.sp <- MDM.mean[-which(rownames(MDM.mean)%in%intersect(rownames(MDM.mean),rownames(MG.mean))),]
tmp.dat <- as.data.frame(apply(MDM.sp,1,function(x){mean(x)}))
colnames(tmp.dat) <- "value"

tmp.dat <- tmp.dat[order(tmp.dat$value),,drop=F]
tmp.dat$pair <- sapply(strsplit(rownames(tmp.dat),"\\."),function(x){x[2]})
tmp.dat$pair <- factor(tmp.dat$pair,levels=as.vector(tmp.dat$pair))
# plot circle barplot 
##################################################
# plot means > 0.5 gene pairs 
#=================================================
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/MDM_naiveCD4.pdf")
ggplot(data = tmp.dat[which(tmp.dat$value>0.5),]) + geom_bar(aes(x=pair,y=value), stat="identity", alpha=0.5)+
	theme_classic()  + scale_fill_viridis_c() +  coord_flip()
dev.off()

































# 2. calculate IFN activatiy in BMDM and MG 
#=====================================================================================================================================================
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM"))
dat <- sub.dat[["RNA"]]@data
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- data.frame(mod)
mod$type <- sub.dat$type.refine

saveRDS(mod,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/LCBM_Myeloid_IFN_mod.RDS")

# Violin Plot 
#=====================================================================================================
library(ggplot2)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/LCBM_Myeloid_IFN_mod.RDS")
dat <- tmp.dat[,c(3,4,7)]

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/IFN_signalingMyeloid.pdf")
ggplot(dat,aes(x=type,y=HALLMARK_INTERFERON_ALPHA_RESPONSE_norm,fill=type))+geom_violin()+ theme_classic()
ggplot(dat,aes(x=type,y=HALLMARK_INTERFERON_GAMMA_RESPONSE_norm,fill=type))+geom_violin()+ theme_classic()
dev.off()









#########################################################################################################
# add a bar plot to show distance of BMDM and T cell
# 2021-1-28
#========================================================================================================

dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CSOmap/BMDM_ALL_CELL_CSOmap.RDS")
coords <- dat[[1]]
cellinfo_tbl <- dat[[2]]
tmp <- aggregate(.~labels,data=cellinfo_tbl[,2:5],FUN=median)

rownames(tmp) <- tmp$labels
tmp$labels<-NULL
as.matrix(dist(tmp)) -> res

# res plot 
#=================================================================
tmp.dist <- res["T_cell",]
tmp.dist <- tmp.dist[-grep("T_cell|unknow",names(tmp.dist))]
tmp.dist <- tmp.dist[order(tmp.dist)]
pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Tcell_distance.pdf")
barplot(tmp.dist,las=2,ylim=c(0,0.15))
dev.off()













###############################################################################################################################
# 2021-3-10
# Treg signature and BMDM MG signature calculate 
#==============================================================================================================================
# GSE14108
dat.BM <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_res.RDS")

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat.BM),c("BMDM_marker","MG_marker","Treg","BMS_test"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- data.frame(mod)
















############################################################################################################################################
# 2021-3-13
# check BMDM and MG IDO1 expression 
#===========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"))
sub.dat@active.ident <- factor(sub.dat$type.refine)

AverageExpression(sub.dat,features=c("IDO1","CD274","IL10","KYNU","QPRT"))



sub.data <- sub.dat[['RNA']]@data
tmp.data <- data.frame(t(as.matrix(sub.data[c("IDO1","CD274","IL10","CXCL10","GBP5"),])))
tmp.data$type <- sub.dat$type.refine
# tmp.data.f <- tmp.data[-which(rowSums(tmp.data[,1:3])==0),]
tmp.res <- aggregate(.~type,data=tmp.data,FUN=mean)
rownames(tmp.res) <- tmp.res$type
tmp.res$type <- NULL

pdf("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Immunosuppress_gene.pdf")
for(i in 1:ncol(tmp.res)){
	barplot(tmp.res[,i],main=colnames(tmp.res)[i],names=c("BMDM","MG"))
}
dev.off()

















library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
TNK <- subset(dat,cells=which(dat$celltype.refine=="Tcell"))

RE_INTE <- function(dat,sample){
    inte.list <- list()
    samplelist <- unique(dat@meta.data[,sample])
    for(i in 1:length(samplelist)){
        tmp <- subset(dat,cells=which(dat@meta.data[,sample]==samplelist[i]))
        DefaultAssay(tmp) <- "RNA"
        inte.list[[i]] <- tmp
    }
    integration.anchors <- FindIntegrationAnchors(object.list = inte.list)
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
    inte <- FindClusters(inte,resolution=1)
    #TSNE
    # if Umap can not use
    inte <- RunTSNE(inte)

    return(inte)
}
TNK <- RE_INTE(TNK,sample="orig.ident")

# saveRDS(TNK,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_TNK.RDS")

library(SingleR)
ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_HumanPrimaryCellAtlasData.RDS")
# ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_NovershternHematopoieticData.RDS")
res.cluster<-SingleR(test=as.matrix(TNK@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=TNK$seurat_clusters,method="cluster")

ref <- readRDS("/public/workspace/lily/software/SingleRData/SingleR_BluePrint.RDS")
res.cluster1<-SingleR(test=as.matrix(TNK@assays$RNA@data),ref=ref,labels=ref$label.fine,clusters=TNK$seurat_clusters,method="cluster")



DefaultAssay(TNK) <- "RNA"
future::plan(multisession, workers=20)
marker <- FindAllMarkers(TNK,assay="RNA",logfc.threshold=0,min.pct=0.05)


TNK$celltype.refine <- "unclassify"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(10,12,15))] <- "NK"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(4,9))] <- "Treg"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(11))] <- "CXCL13T"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(3))] <- "Naive T"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(2,5,6,8))] <- "CD4Tem"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(7,13))] <- "CD8Tcm"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(0,1))] <- "CD8Tcytotoxic"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(17))] <- "CD8Texh"
TNK$celltype.refine[which(TNK$seurat_clusters%in%c(14,16))] <- "proliferationT"

# CD4Tem         CD8Tcm  CD8Tcytotoxic        CD8Texh        CXCL13T
#           6425           2304           4398            431            889
#        Naive T             NK proliferationT           Treg     unclassify



saveRDS(TNK,file="/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_TNK.RDS")


# 1. plot TSNE result 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_TNK.RDS")
cols=c("#8dd593","#b5bbe3","#706357","#023fa5","#f6c4e1","#336600","#f9e498","#774aa4","#cd004b","#c7c7c7")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/DimPlot_TNK.pdf",useDingbats=F)
DimPlot(dat,group.by="celltype.refine",cols=cols,reduction="tsne",raster=F)
dev.off()


# plot supplementary 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/Supplementary_FIg6_FeaturePlot.pdf",useDingbats=F)

FeaturePlot(dat,features=c("CD8A","CD4","CXCL13","MKI67"),label=F,label.size=6,raster=T,reduction="tsne")
FeaturePlot(dat,features=c("CCR7","SELL","TCF7","LEF1"),label=F,label.size=6,raster=T,reduction="tsne") # naive/central memory 
FeaturePlot(dat,features=c("GZMB","GZMK","CTLA4","FOXP3"),label=F,label.size=6,raster=T,reduction="tsne") # effector
FeaturePlot(dat,features=c("HAVCR2","PDCD1","CD3D","TIGIT"),label=F,label.size=6,raster=T,reduction="tsne") # exhausted
dev.off()







# 2. plot celltype percentage in different groups
tmp <- table(dat$type_group,dat$celltype.refine)
tmp.f <- tmp[,-grep("unclassify",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("nMLUAD","MLUAD","LCBM"))
# tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Alveolar Mac","DC","Macrophage","Mast","Monocyte"))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/TNK_celltype_type_group.pdf",useDingbats=F)
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



tmp <- table(dat$type_group,dat$celltype.refine)
rs.list <- list()
for(i in 1:ncol(tmp)){
   rs.list[[i]] <-  c(
    fisher.test(matrix(c(tmp[2,i],(rowSums(tmp)[2]-tmp[2,i]),tmp[2,i],(rowSums(tmp)[3]-tmp[3,i])),ncol=2))$p.value, # MLUAD and nMLUAD
    fisher.test(matrix(c(tmp[1,i],(rowSums(tmp)[1]-tmp[1,i]),tmp[3,i],(rowSums(tmp)[3]-tmp[3,i])),ncol=2))$p.value # LCBM and nMLUAD
    )
    names(rs.list)[i] <- colnames(tmp)[i]
}




# 3. Bulk data verify Treg in LCBM enriched
# cibersort 
# GSE200563 

tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM","MLUNG")),]
dat <- tmp.dat[,rownames(sampleinfo)]
write.table(dat,file="~/software/GSE200563_MLUNG_BM_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)

source("~/software/Cibersort_R.R")
result1 <- CIBERSORT('~/software/LM22.txt','~/software/GSE14108_exp.txt', perm = 100, QN = T)
result2 <- CIBERSORT('~/software/LM22.txt','~/metastasis/data/verify/TCGA_LUAD/LUAD_RNAseq_Exp.txt', perm = 100, QN = T)
result3 <- CIBERSORT('~/software/LM22.txt','~/software/GSE200563_MLUNG_BM_exp.txt', perm = 100, QN = T)

# saveRDS(result2,file="~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_cibersort.RDS")
saveRDS(result1,file="~/metastasis/data/verify/GSE14108/GSE14108_cibsersort.RDS")
saveRDS(result3,file="~/metastasis/data/verify/GSE200563/GSE200563_cibsersort.RDS")

# cibersort result 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/GSE14108_TCGA_treg.cibersort.pdf",useDingbats=F)
boxplot(result1[,9],result2[,9],outline=F,names=c("LCBM","TCGA_LUAD"))
legend("topright",legend=paste0(" pvalue<",0.001))
dev.off()




# 4. use FOXP3 in GSE200563 data
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM","MLUNG") & tmp.sampleinfo$Patient %in%c("Patient_15","Patient_19","Patient_35","Patient_5")),]
dat <- tmp.dat[,rownames(sampleinfo)]
sampleinfo$FOXP3 <- as.numeric(dat["FOXP3",])

sampleinfo$group <- factor(sampleinfo$group,levels=c("MLUNG","BM"))
p <- ggplot(sampleinfo, aes(x = group, y = FOXP3)) +
geom_boxplot(aes(fill = group), show.legend = FALSE, width = 0.6) +  #绘制箱线图
scale_fill_manual(values = c('#FE7280', '#AC88FF')) +  #箱线图的填充色
geom_point(size = 2) +  #绘制样本点
geom_line(aes(group = Patient), color = 'black', lwd = 0.5) +  #绘制配对样本间连线
##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black', size = 1), panel.background = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
    axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 20, color = 'black')) +
labs(x = 'Group', y = 'expression of FOXP3', title = 'FOXP3', subtitle = 'MLUNG vs. BM')

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/GSE200563_FOXP3.paired.pdf",useDingbats=F)
print(p)
dev.off()








#==================================================================================================================================================
# use cellphoneDB result 
# 1. calculate communication numbers 
sampledir <- dir("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/")
rs.list <- list()
for(i in 1:length(sampledir)){
    tmp.dat <- read.table(paste0("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/",sampledir[i],"/significant_means.txt"),sep="\t",header=T)
    rownames(tmp.dat) <- tmp.dat$id_cp_interaction

    tmp.f <- tmp.dat[,c("Tumor.Treg","Macrophage_P.Treg","Macrophage_N.Treg","MG.Treg")]
    tmp.res <- tmp.f[which(rownames(tmp.f)%in%names(which(apply(tmp.f,1,function(x){if(all(is.na(x))){"T"}else{"F"}})=="F"))),]

    tmp.num <- apply(tmp.res,2,function(x){length(which(!is.na(x)))})

    rs.list[[i]] <- tmp.num
    names(rs.list)[i] <- sampledir[i]

}

rs.mat <- matrix(unlist(rs.list),ncol=4,byrow=T)
colnames(rs.mat) <- c("Tumor.Treg","Macrophage_P.Treg","Macrophage_N.Treg","MG.Treg")
rownames(rs.mat) <- names(rs.list)

plotdat <- as.data.frame(apply(rs.mat,2,function(x){mean(x)}))
plotdat$group <- rownames(plotdat)
plotdat$se <- (apply(rs.mat,2,function(x){mean(x)}))/sqrt(6)
colnames(plotdat)[1] <- "mean"

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/scRNA_LCBM_CellphoneDB_Treg_num.pdf",useDingbats=F)
ggplot(plotdat) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-se, ymax=mean+se), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme_bw()
dev.off()








#=======================================================================================================================
# 2023-02-20
sampledir <- dir("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/")
rs.list <- list()
for(i in 1:length(sampledir)){
    tmp.dat <- read.table(paste0("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/",sampledir[i],"/significant_means.txt"),sep="\t",header=T)
    rownames(tmp.dat) <- tmp.dat$id_cp_interaction

    tmp.f <- tmp.dat[,grep("Tumor$|unclassify$|MG$|Macrophage_P$|NK$|Macrophage_N$",grep("^Tumor.*",colnames(tmp.dat),value=T),value=T,invert=T)]
    tmp.res <- tmp.f[which(rownames(tmp.f)%in%names(which(apply(tmp.f,1,function(x){if(all(is.na(x))){"T"}else{"F"}})=="F"))),]

    tmp.num <- apply(tmp.res,2,function(x){length(which(!is.na(x)))})

    rs.list[[i]] <- tmp.num
    names(rs.list)[i] <- sampledir[i]

}

median(sapply(rs.list,function(x){sum(x)}))

# P 86       329       284       293       226       164
# 255
# N 59       216       155       153       114       102
# 133.5
# Tumor 47       282       233       206       135       154
# 180
# MG 60       172       190       267       131       134
# 153 
# interaction numbers to all T cells
IntNum <- data.frame(MacP=c(86,329,284,293,226,164),MacN=c(59,216,155,153,114,102),Tumor=c(47,282,233,206,135,154),MG=c(60,172,190,267,131,134),stringsAsFactors=F)
rownames(IntNum) <- c("A20190305","A20190312","D0927","E0927","Pair_BM","T_Bsc1")

plotdat <- reshape2::melt(IntNum)
colnames(plotdat) <- c("Group","Num")
# no group significant
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/scRNA_LCBM_CellphoneDB_Tcell_num.pdf",useDingbats=F)
ggplot(plotdat,aes(x=Group,y=Num,fill=Group))+
	geom_bar(fun="median",position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = plotdat, aes(x=Group,y = Num,fill=Group),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c("#d4c78c","#5ec6f2","#003468","#4b1702"))+ theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
dev.off()
# boxplot(value~variable,data=plotdat)
# beeswarm::beeswarm(value~variable,data=plotdat,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)









#========================================================================================================================================
# which max for Macrophage_P and T cells
sampledir <- dir("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/")
rs.list <- list()
for(i in 1:length(sampledir)){
    tmp.dat <- read.table(paste0("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/",sampledir[i],"/significant_means.txt"),sep="\t",header=T)
    rownames(tmp.dat) <- tmp.dat$id_cp_interaction

    tmp.f <- tmp.dat[,grep("Tumor$|unclassify$|MG$|Macrophage_P$|NK$|Macrophage_N$",grep("^Macrophage_P.*",colnames(tmp.dat),value=T),value=T,invert=T)]
    tmp.res <- tmp.f[which(rownames(tmp.f)%in%names(which(apply(tmp.f,1,function(x){if(all(is.na(x))){"T"}else{"F"}})=="F"))),]

    tmp.num <- apply(tmp.res,2,function(x){length(which(!is.na(x)))})

    rs.list[[i]] <- tmp.num
    names(rs.list)[i] <- sampledir[i]

}
rs.list$Pair_BM["Macrophage_P.CD8Texh"] <- 0
rs.list$Pair_BM <- rs.list$Pair_BM[names(rs.list$E0927)]
rs.mat <- matrix(unlist(rs.list),ncol=8,byrow=T)
colnames(rs.mat) <- names(rs.list$T_Bsc1)
rownames(rs.mat) <- names(rs.list)


tmpdat <- matrix(0,ncol=9,nrow=9)
tmpdat[1,2:9] <- apply(rs.mat,2,function(x){mean(x)})
colnames(tmpdat)[2:9] <- sapply(strsplit(colnames(rs.mat),"\\."),function(x){x[[2]]}) 
rownames(tmpdat)[2:9] <- sapply(strsplit(colnames(rs.mat),"\\."),function(x){x[[2]]}) 
rownames(tmpdat)[1] <- "Macrophage_P"
colnames(tmpdat)[1] <- "Macrophage_P"

library(CellChat)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/scRNA_LCBM_CellphoneDB_MacrophageP_Tcell_num.pdf",useDingbats=F)
netVisual_circle(net=tmpdat,weight.scale=T,label.edge=T,edge.width.max=20)
dev.off()















#================================================================================================================
# specific interaction pair
# 2023-02-20
sampledir <- dir("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/")
rs.list <- list()
for(i in 1:length(sampledir)){
    tmp.dat <- read.table(paste0("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/",sampledir[i],"/significant_means.txt"),sep="\t",header=T)
    rownames(tmp.dat) <- tmp.dat$id_cp_interaction

    rs.list[[i]] <- as.character(unlist(tmp.dat[which(!is.na(tmp.dat[,"Macrophage_P.proliferationT"])),2]))
    names(rs.list)[i] <- sampledir[i]

}

a <- names(which(table(unlist(rs.list))==6))

exp.list <- list()
for(i in 1:length(sampledir)){
    tmp.dat <- read.table(paste0("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/",sampledir[i],"/significant_means.txt"),sep="\t",header=T)
    rownames(tmp.dat) <- tmp.dat$id_cp_interaction
    tmp.exp <- tmp.dat[which(tmp.dat[,2]%in%a),"Macrophage_P.proliferationT"]
    names(tmp.exp) <- tmp.dat[which(tmp.dat[,2]%in%a),2]
    exp.list[[i]] <- tmp.exp
    names(exp.list)[i] <- sampledir[i]

}

exp.mat <- matrix(unlist(exp.list),ncol=8,byrow=T)
pairs <- apply(exp.mat,2,function(x){mean(x)})
names(pairs) <- a

pdat <- data.frame(source=sapply(strsplit(names(pairs),"_"),function(x){x[[1]]}),
    target=sapply(strsplit(names(pairs),"_"),function(x){x[[2]]}),
    weight=pairs)

nodes <- data.frame(name=c(as.character(pdat$source), as.character(pdat$target)) %>% unique())
pdat$IDsource <- match(pdat$source, nodes$name)-1
pdat$IDtarget <- match(pdat$target, nodes$name)-1
save(pdat,nodes,file="/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/scRNA_LCBM_CellphoneDB_MacrophageP_proliferationT.RData")


library(networkD3)
sankeyNetwork(Links = pdat, Nodes = nodes, Source = "IDsource", Target = "IDtarget",
Value = "weight", NodeID = "name",nodeWidth =10,units = 'TWh',
height=300,width=300,colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),fontSize = 8)

webshot::webshot("C:\\Users\\Lily\\Oned_lly\\OneDrive - njmu.edu.cn\\version_12_23\\Download_plot\\Fig6/Sankeyplot_res.html",
    "C:\\Users\\Lily\\Oned_lly\\OneDrive - njmu.edu.cn\\version_12_23\\Download_plot\\Fig6/sankey.pdf")





##########################################################################################################################################
# 2023-02-20
# use GSE200563 to verify

tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM")),]
dat <- tmp.dat[,rownames(sampleinfo)]

immune <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_cibsersort.RDS")

# 不好看
# cor.test(as.numeric(dat["PTGER4",]),as.numeric(immune[,9]),method="spearman")

# plot(as.numeric(dat["PTGER4",]),as.numeric(immune[,9]),main=" GSE200563 ")
# abline(lm(as.numeric(immune[,9])~as.numeric(dat["PTGER4",])),col="red")
# legend("topright",legend=paste0("rho= 0.28"," P = 0.046","N=51"))

# dev.off()


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/scRNA_LCBM_CellphoneDB_PTGER4_PDCD1_cor.pdf",useDingbats=F)

cor.test(as.numeric(dat["PTGER4",]),as.numeric(dat["PDCD1",]),method="spearman")
plot(as.numeric(dat["PTGER4",]),as.numeric(dat["PDCD1",]),main=" GSE200563 ")
abline(lm(as.numeric(dat["PDCD1",])~as.numeric(dat["PTGER4",])),col="red")
legend("topright",legend=paste0("rho= 0.36"," P = 0.01","N=51"))

dev.off()



# # PTGER4
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("Tcell_dys","Tcell_act"),"/public/workspace/lily/Lung2Brain/Version6/Caculate_Data/MOD_file/",permN=0)
# mod <- as.data.frame(mod)
# mod$score <- mod$Tcell_dys_norm- mod$Tcell_act_norm
# mod$PTGER4 <- as.numeric(dat["PTGER4",])




# immune checkpoint signature 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4334715/
# PD-1, LAG-3, TIM-3, VISTA, CD244, CD160, and BTLA
# genes <- c("PDCD1","LAG3","TIM3","VISTA","CD244","CD160","BTLA")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(genes,"immune_checkpoint",out="/public/workspace/lily/MOD_file/immune_checkpoint.mod") # make a mod file 

source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(dat,c("immune_checkpoint"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$PTGER4 <- as.numeric(dat["PTGER4",])
cor.test(mod$PTGER4,mod$immune_checkpoint_norm,method="spearman")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig6/GSE200563_BM_PTGER4_immuneCheckpoint_cor.pdf",useDingbats=F)

cor.test(mod$PTGER4,mod$immune_checkpoint_norm,method="spearman")
plot(mod$PTGER4,mod$immune_checkpoint_norm,main=" GSE200563 ")
abline(lm(mod$immune_checkpoint_norm~mod$PTGER4),col="red")
legend("topright",legend=paste0("rho= 0.53"," P = 0.007","N=24"))

dev.off()






# 2021-12-30
# this program is used to analysis GSE131907 data 
#=========================================================================================================================================
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")

# just get satge I and stage II
# get 9 samples
# early <- subset(tmp.dat,cells=which(tmp.dat$Sample%in%c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19",
#     "LUNG_T20","LUNG_T25",'LUNG_T30','LUNG_T34')))

for(i in c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T25",'LUNG_T30','LUNG_T34')){
    tmp_dat <- subset(tmp.dat,cells=which(tmp.dat$Sample==i))
    # normalize and filter 
# seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	tmp_dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp_dat, pattern = "^MT-")
	# pdf(paste0("/public/workspace/lily/Lung2Brain/GSE131907/early/Plot/",i, "_vlnPlot_prepare.pdf"))
	# p<-VlnPlot(object = tmp_dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
	# print(p)
    # dev.off()
    saveRDS(tmp_dat,file=paste0("/public/workspace/lily/Lung2Brain/GSE131907/early/",i, ".RDS"))
}



















# 2021-12-30 
# before integration ,try to find signature 
#============================================================================================================================================
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/GSE131907/all_cell.RDS")

# just get satge I and stage II
# get 9 samples
early <- subset(tmp.dat,cells=which(tmp.dat$Sample%in%c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19",
    "LUNG_T20","LUNG_T25",'LUNG_T30','LUNG_T34')))

	tmp_dat = NormalizeData(object =early)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	early <- ScaleData(object = tmp_dat, features = all.genes)

tumor <- subset(early,cells=which(early$Cell_subtype%in%c("tS1","tS2","tS3")))

# re inte 
inte.list <- list()
samplelist <- unique(tumor$Sample)
for(i in 1: length(samplelist)){
    tmp <- subset(tumor,cells=which(tumor$Sample==samplelist[i]))
    DefaultAssay(tmp) <- "RNA"
    inte.list[[i]] <- tmp
}
integration.anchors <- FindIntegrationAnchors(object.list = inte.list,k.filter=50)
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
inte <- FindClusters(inte,resolution=0.5)
#TSNE
# if Umap can not use
inte <- RunTSNE(inte)
inte <- RunUMAP(inte,dims=1:10)


saveRDS(tumor,file="/public/workspace/lily/Lung2Brain/Version6/Data/SanSun_early_tumor.RDS")







































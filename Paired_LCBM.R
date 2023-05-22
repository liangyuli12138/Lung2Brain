
# 2021-12-15
# use this program to analysis paired LCBM data 
# this data run in .3 
# and RDS object will trans to .5 
#===============================================================================================================================================
library(Seurat)

# LCBM data 
tmp <- Read10X(data.dir = "/public/workspace/lily/LCBM_pair/R21125541/outs/filtered_feature_bc_matrix/")
tmp.dat <- CreateSeuratObject(counts = tmp,  project = "Pair_BM",min.cells = 3, min.features = 200)

prepare <- function(tmp_dat){
    tmp_dat[["percent.mt"]] <- PercentageFeatureSet(object = tmp_dat, pattern = "^MT-")
    tmp_dat = subset(x=tmp_dat,subset=nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
# seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
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
	return(tmp_dat)
}

dat <- prepare(tmp.dat)

saveRDS(dat,file=)




# Lung cancer data
tmp <- Read10X(data.dir = "/public/workspace/lily/LCBM_pair/R21136163/outs/filtered_feature_bc_matrix/")
tmp.dat <- CreateSeuratObject(counts = tmp,  project = "Pair_Lung",min.cells = 3, min.features = 200)
dat <- prepare(tmp.dat)


































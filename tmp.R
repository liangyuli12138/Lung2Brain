
samplefile <- c("GSM3516664_MSK_LX666_METASTASIS_dense.csv.gz","GSM3516668_MSK_LX255B_METASTASIS_dense.csv.gz","GSM3516671_MSK_LX681_METASTASIS_dense.csv.gz","GSM3516677_MSK_LX699_METASTASIS_dense.csv.gz","GSM3516678_MSK_LX701_METASTASIS_dense.csv.gz")
dat.list <- list()
for(i in 1:length(samplefile)){
    tmp.dat <- read.csv(samplefile[i])
    samplename <- paste0("S",substr(strsplit(samplefile[i],"_")[[1]][1],9,10))
    rownames(tmp.dat) <- paste0(samplename,tmp.dat$X)
    tmp.dat$X <- NULL

    tmp_dat<- CreateSeuratObject(counts = t(tmp.dat),  project = samplename,min.cells = 3, min.features = 200)
    # seurat object
	tmp_dat = NormalizeData(object = tmp_dat)
	tmp_dat <- FindVariableFeatures(object = tmp_dat)
	# scaling
	all.genes <- rownames(x = tmp_dat)
	tmp_dat <- ScaleData(object = tmp_dat, features = all.genes)
    tmp_dat <- RunPCA(object = tmp_dat, features = VariableFeatures(object = tmp_dat))
	# clustering
	tmp_dat <- FindNeighbors(object = tmp_dat,dims=1:10)
	# select proper resolution
	tmp_dat <- FindClusters(object = tmp_dat,resolution=0.5)
    tmp_dat <- RunTSNE(object = tmp_dat,dims=1:10,check_duplicates = FALSE)
	tmp_dat <- RunUMAP(tmp_dat,dims=1:10)
    dat.list[[i]] <- tmp_dat
}

pdf("~/tmp.pdf")
for(i in 1:length(dat.list)){
    p <- FeaturePlot(dat.list[[i]],features=c("EGFR","EPCAM","OSMR","OSM"),label=T,label.size=8)
    tmp.dat <- dat.list[[i]]
    tmp.dat$OSMR <- ifelse(tmp.dat[["RNA"]]@data["OSMR",]>0,"Y","N")
    print(apply(table(tmp.dat$OSMR,tmp.dat$seurat_clusters),1,function(x){x/sum(x)}))
    print(mean(tmp.dat[["RNA"]]@data["OSMR",]))
    cat("\n")
    print(p)
}
dev.off()



































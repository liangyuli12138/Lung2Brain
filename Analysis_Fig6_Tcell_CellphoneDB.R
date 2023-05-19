
# 2022-8-25
# do cellphoneDB for Macrophage Tumor and T cells in LCBM 
#===============================================================================================================================
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
tumor <- subset(tmp.dat,cells=which(tmp.dat$celltype.refine=="Tumor"&tmp.dat$type_group=="LCBM"))
tumor$CCC.type <- "Tumor"

tmp.mye <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_Myeloid.RDS")
mac <- subset(tmp.mye,cells=which(tmp.mye$celltype.refine%in%c("Macrophage","MG") & tmp.mye$type_group=="LCBM"))
mac$OSM <- ifelse(mac[["RNA"]]@data["OSM",]>0,"P","N")
mac$CCC.type <- paste0(mac$celltype.refine,"_",mac$OSM)
mac$CCC.type[which(mac$celltype.refine =="MG")] <- "MG"

tmp.t <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16_TNK.RDS")
T <- subset(tmp.t,cells=which(tmp.t$type_group=="LCBM"))
T$CCC.type <- T$celltype.refine

dat <- merge(x=tumor,y=c(mac,T))

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

for(i in unique(dat$orig.ident)){
    subdat <- subset(dat,cells=which(dat$orig.ident==i))
    dir.create(paste0("/public/workspace/lily/Lung2Brain/Version6/Data/CellphoneDB/",i))
    outpath=paste0("/public/workspace/lily/Lung2Brain/Version6/Data/CellphoneDB/",i)
    cellphoneDB_input(subdat,"CCC.type",paste0(outpath,"/expr.txt"),paste0(outpath,"/cellinfo.txt"))

}




# run in shell
# for i in `ls /public/workspace/lily/Lung2Brain/Version6/Data/CellphoneDB/`
# do
#     echo $i
# 	bytlib load lib64.pool
# 	bytlib load sqlite3-snapshot
# 	bytlib load python-3.6.6
# 	bytlib load cellphonedb-2.1.1

# 	mkdir -p /public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/${i}/
# 	RESPATH="/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/${i}/"
# 	# run CellphoneDB
# 	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/Version6/Data/CellphoneDB/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/Version6/Data/CellphoneDB/${i}/expr.txt --threads 8 \
# 	--output-path=${RESPATH} --counts-data gene_name

# done





#============================================================================================================================
# 2022-8-26
# analysis with result 

tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/A20190305/significant_means.txt",sep="\t",header=T)
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/A20190312/significant_means.txt",sep="\t",header=T)

tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/T_Bsc1/significant_means.txt",sep="\t",header=T)
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/Pair_BM/significant_means.txt",sep="\t",header=T)
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/E0927/significant_means.txt",sep="\t",header=T)
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Version6/CellphoneDB/Res/D0927/significant_means.txt",sep="\t",header=T)

tmp.dat[which(!is.na(tmp.dat[,"Tumor.Treg"])),2]
tmp.dat[which(!is.na(tmp.dat[,"Macrophage_P.proliferationT"])),2]
tmp.dat[which(!is.na(tmp.dat[,"Macrophage_N.Treg"])),2]
tmp.dat[which(!is.na(tmp.dat[,"MG.Treg"])),2]





































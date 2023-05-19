
#========================================================================================
# 2020-11-30 
# should use more pure Myeloid cell to run CellphoneDB 
#========================================================================================
library(Seurat)
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_MG_mod.RDS")
# use pvalue to filter 
tmp.type <- as.data.frame(apply(mod[,5:6],1,function(x){
    tmp <- length(names(which(x[]<0.05)))
	if(tmp==0|tmp>1){
		c("NA")
	}else{
		names(which(x[]<0.05))
	}
}))
colnames(tmp.type) <- "mod_res"
Myeloid$type.CCC <- tmp.type$mod_res
Myeloid$type.CCC <- gsub("_marker_pval","",as.vector(Myeloid$type.CCC))
saveRDS(Myeloid,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
#====================================
# table(Myeloid$type_group,Myeloid$type.CCC) -> tmp
# tmp.f <- tmp[,c(1,2)]


#========================================================================================
# ready to run CellphoneDB 
#========================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type%in%c("malignant","Myeloid")))
sub.dat$type.CCC <- sub.dat$type 
sub.dat$type.CCC[which(colnames(sub.dat)%in%colnames(Myeloid))] <- Myeloid$type.CCC
sub.dat.f <- subset(sub.dat,cells=which(sub.dat$type.CCC%in%c("BMDM","malignant","MG")))
sub.dat.f$sample <- sub.dat.f$orig.ident
sub.dat.f$sample[grep("RD-20180817-001-SR18271",sub.dat.f$sample)] <- "lesion1"
sub.dat.f$sample[grep("RD-20180817-002-SR18271",sub.dat.f$sample)] <- "lesion2"
sub.dat.f$sample[grep("T-Bsc1",sub.dat.f$sample)] <- "T_Bsc1"

#===========================
# make output 
#===========================
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

sample_name <- unique(sub.dat.f$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/",sample_name[i],"/")
    tmp <- subset(sub.dat.f,cells=which(sub.dat.f$sample==sample_name[i]))
    cellphoneDB_input(tmp,"type.CCC",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

#============================================
# run in shell
#============================================

for i in `ls /public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done

#========================================================================================================================
# run result 
# pvalue 
#========================================================================================================================
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/")
Malignant.BMDM <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","malignant.BMDM")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".malignant.BMDM") 
	Malignant.BMDM[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Malignant.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Malignant.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) # get co_id for intersection pair id 

Malignant.BMDM.res <- cbind(Malignant.BMDM[[1]][which(Malignant.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Malignant.BMDM[[2]][which(Malignant.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Malignant.BMDM)){
	Malignant.BMDM.res <- cbind(Malignant.BMDM.res,Malignant.BMDM[[i]][which(Malignant.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.BMDM.res <- Malignant.BMDM.res[,c(1,grep("BMDM",colnames(Malignant.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.BMDM.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair
#=======================================================================================================================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]
Com.res.f <- Com.res.f[,c("A20190305.malignant.BMDM","A20190312.malignant.BMDM","T_Bsc1.malignant.BMDM","lesion1.malignant.BMDM","lesion2.malignant.BMDM",
	"BT1296.malignant.BMDM","BT1297.malignant.BMDM","scrBT1431m.malignant.BMDM","scrBT1432m.malignant.BMDM")]
# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Malignant.BMDM.pvalue.f.RDS")
















#=====================================================================================================================
# means 
#=====================================================================================================================
rm(list=ls())
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/")
Malignant.BMDM <- list()
for(i in 1:length(sample_name)){
	#tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.mean.f <- tmp.mean[,c("id_cp_interaction","malignant.BMDM")]
	colnames(tmp.mean.f)[2] <- paste0(sample_name[i],".malignant.BMDM") 
	Malignant.BMDM[[i]] <- tmp.mean.f
}
#===========================================================================
tmp.id <- c()
for(i in 1:length(Malignant.BMDM)){
	tmp.id <- c(tmp.id,as.vector(Malignant.BMDM[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==9)) 
#===============================================
Malignant.BMDM.res <- cbind(Malignant.BMDM[[1]][which(Malignant.BMDM[[1]]$id_cp_interaction%in%co_id),],
							Malignant.BMDM[[2]][which(Malignant.BMDM[[2]]$id_cp_interaction%in%co_id),])
# get co_id means matrix 
for(i in 3:length(Malignant.BMDM)){
	Malignant.BMDM.res <- cbind(Malignant.BMDM.res,Malignant.BMDM[[i]][which(Malignant.BMDM[[i]]$id_cp_interaction%in%co_id),])
}
# do some modify
Malignant.BMDM.res <- Malignant.BMDM.res[,c(1,grep("BMDM",colnames(Malignant.BMDM.res)))]
# add pair gene information 
Communcation.res <- merge(Malignant.BMDM.res,tmp.mean[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair
Com.res <- Communcation.res[,grep("BMDM",colnames(Communcation.res))]
Com.res.f <- Com.res[,c("A20190305.malignant.BMDM","A20190312.malignant.BMDM","T_Bsc1.malignant.BMDM","lesion1.malignant.BMDM","lesion2.malignant.BMDM",
	"BT1296.malignant.BMDM","BT1297.malignant.BMDM","scrBT1431m.malignant.BMDM","scrBT1432m.malignant.BMDM")]
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Malignant.BMDM.mean.f.RDS")


#===============================
pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Malignant.BMDM.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/Malignant.BMDM.mean.f.RDS")
means.res <- means.res[-which(rowSums(means.res)<9),]

# should make a id to merge different samples 
library(reshape)
pvalue <- melt(as.matrix(pvalue.res),id="col.names")
colnames(pvalue) <- c("gene_pair","samples","pvalue")
pvalue$id <- paste0(pvalue$gene_pair,".",pvalue$samples)
# should make a id to merge different samples
mean <- melt(as.matrix(means.res),id="col.names")
colnames(mean) <- c("gene_pair","samples","means_value")
mean$id <- paste0(mean$gene_pair,".",mean$samples)
#========================================================================
merge(mean,pvalue,by="id") -> res.final
res.final$Samples <- sapply(as.vector(res.final$samples.x),function(x){strsplit(x,"\\.")[[1]][1]})
res.final$Samples <- factor(res.final$Samples,levels=c("lesion1","lesion2","A20190305","A20190312","T_Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))

res.final.f <- res.final
res.final.f$means_value[which(res.final.f$pvalue==0)] <- 0


pdf("/public/workspace/lily/Lung2Brain/TME/result_plot/Malignant.BMDM_bubble.pdf",useDingbats=F)

library(ggplot2)
ggplot(res.final.f, aes(x=Samples, y=gene_pair.y, size=pvalue)) + geom_point(aes(colour = means_value))+
	theme_bw()+theme(axis.text.x = element_text(angle=45,size=7,hjust=0.4))+
	#scale_color_distiller(palette = "Greens")+
	#scale_color_gradient(low = "white",high = "ForestGreen")+
	#scale_colour_gradientn(low = "ForestGreen",mid="white",high = "#F29403")+
	#scale_colour_gradientn(colors=c('white','#F29403'),values=c(0,0.5,20))+
	scale_colour_gradientn(colours = colorRampPalette(c('white','#F29403'))(500))+
	theme(panel.grid.major = element_blank())

dev.off()



























#========================================
# multiple lesion test 
# 2020-11-3
#========================================

file <- dir("/public/workspace/lily/Lung2Brain/multiple_lesion_res")[grep("count",dir("/public/workspace/lily/Lung2Brain/multiple_lesion_res"))]
count.list <- list()
for(i in 1: length(file)){
    tmp <- read.table(paste0("/public/workspace/lily/Lung2Brain/multiple_lesion_res/",file[i]))
    colnames(tmp) <- c("gene_ensemble",file[i])
    assign(file[i],tmp)
    count.list[[i]] <- tmp
}

res <- Reduce(function(x,y){merge(x,y,by="gene_ensemble")},count.list,accumulate =FALSE)

res.f <- res[grep("ENSG",res$gene_ensemble),]

#========================================
# load ensemble info 
#========================================
dat <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/hg19_ensemble.txt",header=T,sep="\t")
dat$length <- dat$Gene.end..bp. -dat$Gene.start..bp.
colnames(dat)[1] <- "gene_ensemble"

res.tmp <- merge(res.f,dat,by="gene_ensemble")
rownames(res.tmp) <- res.tmp$gene_ensemble
res.tmp[,c(2:9,13)] -> mat
mat.tmp <- apply(mat,1,function(x){x/x[length(x)]})
mat.res <- apply(mat.tmp,1,function(x){(x/sum(x))*10^6})
mat.res <- data.frame(mat.res)
mat.res$gene_ensemble <- rownames(mat.res)
res <- merge(mat.res,dat[,c(1,2)],by="gene_ensemble")
res$gene_ensemble <- NULL
res$length <- NULL

aggregate(.~Gene.name,data=res,FUN=median) -> res.final
rownames(res.final) <- res.final$Gene.name
res.final$Gene.name <- NULL
head(res.final)

saveRDS(res.final,file="/public/workspace/lily/Lung2Brain/multiple_lesion_res/Expr.RDS")






#=====================================================================================================
# calculate ssgsea between multiple lesion 
#=====================================================================================================
dat <- readRDS("/public/workspace/lily/Lung2Brain/multiple_lesion_res/Expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)

#=====================================================================================================
# calculate mod is not ok, there are some big difference between different lesion 
#
#=====================================================================================================
library(pheatmap)
gene <- readRDS("/public/workspace/lily/Lung2Brain/inte7/BMS_gene.RDS")
pheatmap(dat[gene,],scale="row")































# 2021-10-13
# also calculate three groups gene expression 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
dat@active.ident <- factor(dat$group)

pdf("/public/workspace/lily/Lung2Brain/Version5/Brain_like/LCBM_gene_number.pdf",useDingbats=F)
boxplot(nFeature_RNA~group,data=dat@meta.data,FUN=median,main="expressed gene numbers",outline=F)
dev.off()


# GSE123902 data 
rm(list=ls())
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig5/GSE123902_BM_tumor.RDS")
pdf("/public/workspace/lily/Lung2Brain/Version5/Brain_like/GSE123902_gene_number.pdf",useDingbats=F)
boxplot(nFeature_RNA~tmp.group,data=dat@meta.data,FUN=median,main="expressed gene numbers",outline=F)
dev.off()




# tumor.d and tumor.e  data 
rm(list=ls())
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.E.RDS")
pdf("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Tumor_E_gene_number.pdf",useDingbats=F)
boxplot(nFeature_RNA~tmp.group,data=dat@meta.data,FUN=median,main="expressed gene numbers",outline=F)
dev.off()

rm(list=ls())
dat <- readRDS("/public/workspace/lily/Lung2Brain/Multiple_LungBrain/tumor.D.RDS")
pdf("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Tumor_D_gene_number.pdf",useDingbats=F)
boxplot(nFeature_RNA~tmp.group,data=dat@meta.data,FUN=median,main="expressed gene numbers",outline=F)
dev.off()






#============================================================================================================================================================
# calculate gene permutation for LCBM subtumor 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")

# calculate sub group gene correlation 
lung <- subset(dat,cells=which(dat$group=="group09"))
lung.dat <- lung[["RNA"]]@data
tmp.percent <- apply(as.matrix(lung.dat),1,function(x){length(which(x>0))/ncol(lung.dat)})
lung.dat.f <- lung.dat[-which(tmp.percent<0.05),] # filter a lot of gene 
# then calculate 
library(Hmisc)
tmp.res <- rcorr(t(as.matrix(lung.dat.f)))
res.f <- tmp.res$r
res.f[!upper.tri(res.f, diag = FALSE)]=0 # tri matrix 
saveRDS(res.f,file="/public/workspace/lily/Lung2Brain/Version5/Brain_like/Lung.cor.mat.RDS")
#======================================================================================================
brain <- subset(dat,cells=which(dat$group=="group16"))
brain.dat <- brain[["RNA"]]@data
tmp.percent <- apply(as.matrix(brain.dat),1,function(x){length(which(x>0))/ncol(brain.dat)})
brain.dat.f <- brain.dat[-which(tmp.percent<0.05),] # filter a lot of gene 
# then calculate 
library(Hmisc)
tmp.res <- rcorr(t(as.matrix(brain.dat.f)))
res.f <- tmp.res$r
res.f[!upper.tri(res.f, diag = FALSE)]=0 # tri matrix 
saveRDS(res.f,file="/public/workspace/lily/Lung2Brain/Version5/Brain_like/Brain.cor.mat.RDS")
#======================================================================================================
stem <- subset(dat,cells=which(dat$group=="group47"))
stem.dat <- stem[["RNA"]]@data
tmp.percent <- apply(as.matrix(stem.dat),1,function(x){length(which(x>0))/ncol(stem.dat)})
stem.dat.f <- stem.dat[-which(tmp.percent<0.05),] # filter a lot of gene 
# then calculate 
library(Hmisc)
tmp.res <- rcorr(t(as.matrix(stem.dat.f)))
res.f <- tmp.res$r
res.f[!upper.tri(res.f, diag = FALSE)]=0 # tri matrix 
saveRDS(res.f,file="/public/workspace/lily/Lung2Brain/Version5/Brain_like/Stem.cor.mat.RDS")




#============================================================================================================================================
# now random select cells to construct net work 
#============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Fig4/Data/LCBM_sub.tumor.RDS")
stem <- subset(dat,cells=which(dat$group=="group47"))










#=============================================================================================================================================
# get gene and construct net 
# lung 
mat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Lung.cor.mat.RDS")
lung.res <- data.frame(gene1="A",gene2="B",correlation=1,stringsAsFactors=F)
threshold <- 0.3
for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
        if(mat[i,j]>threshold){
            lung.res <- rbind(lung.res,c(rownames(mat)[i],colnames(mat)[j],mat[i,j]))
        }
        
    }
}

# set network 
net.dat <- lung.res[-1,]
colnames(net.dat) <- c("from","to","correlation")
net <- graph_from_data_frame(net.dat, directed=TRUE)
write.table(net.dat,file=paste0("/public/workspace/lily/Lung2Brain/Version5/Brain_like/","Lung_cor_",threshold,"_table.txt"),sep="\t",row.names=F,col.names=T,quote=F)

library(igraph)
# degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)






# brain 
mat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Brain.cor.mat.RDS")
#length(grep("^RPL|^RPS",colnames(mat)))/ncol(mat)
brain.res <- data.frame(gene1="A",gene2="B",correlation=1,stringsAsFactors=F)
threshold <- 0.3
for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
        if(mat[i,j]>threshold){
            brain.res <- rbind(brain.res,c(rownames(mat)[i],colnames(mat)[j],mat[i,j]))
        }
        
    }
}

# set network 
library(igraph)
net.dat <- brain.res[-1,]
dim(net.dat)
colnames(net.dat) <- c("from","to","correlation")
net <- graph_from_data_frame(net.dat, directed=TRUE)
write.table(net.dat,file=paste0("/public/workspace/lily/Lung2Brain/Version5/Brain_like/","Brain_cor_",threshold,"_table.txt"),sep="\t",row.names=F,col.names=T,quote=F)


# degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)






# stem like 
mat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Stem.cor.mat.RDS")
stem.res <- data.frame(gene1="A",gene2="B",correlation=1,stringsAsFactors=F)
threshold <- 0.3
for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
        if(mat[i,j]>threshold){
            stem.res <- rbind(stem.res,c(rownames(mat)[i],colnames(mat)[j],mat[i,j]))
        }
        
    }
}

# set network 
library(igraph)
net.dat <- stem.res[-1,]
dim(net.dat)
colnames(net.dat) <- c("from","to","correlation")
net <- graph_from_data_frame(net.dat, directed=TRUE)
write.table(net.dat,file=paste0("/public/workspace/lily/Lung2Brain/Version5/Brain_like/","Stem_cor_",threshold,"_table.txt"),sep="\t",row.names=F,col.names=T,quote=F)


# degree(net, mode="in")
centr_degree(net, mode="in", normalized=T)
transitivity(net, type = "average")















#===========================================================================================================================================
# random gene numbers to construct net work 
# Lung like tumor cells
lung.mat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/Brain_like/Lung.cor.mat.RDS")
idx <- sample(1:ncol(lung.mat),size=4315,replace=F)
mat <- lung.mat[idx,idx]

# calculate net work 
lung.res <- data.frame(gene1="A",gene2="B",correlation=1,stringsAsFactors=F)
threshold <- 0.5
for(i in 1:nrow(mat)){
    for(j in 1:ncol(mat)){
        if(mat[i,j]>threshold){
            lung.res <- rbind(lung.res,c(rownames(mat)[i],colnames(mat)[j],mat[i,j]))
        }
        
    }
}

# set network 
library(igraph)
net.dat <- lung.res[-1,]
colnames(net.dat) <- c("from","to","correlation")
net <- graph_from_data_frame(net.dat, directed=TRUE)
transitivity(net, type = "average")



















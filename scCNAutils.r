
#========================================================================================
# 2020-12-2 
# infercnv use a new way 
#========================================================================================
library(scCNAutils)
library(Seurat)

# just a test 
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
tmp <- subset(dat,cells=which(dat$orig.ident=="A20190312"))
sub.dat <- subset(tmp,cells=which(tmp$type%in%c("B_cell","maliganant")))
Bcell <- colnames(sub.dat)[which(sub.dat$type=="B_cell")]

#========================================================================================
# run 

expr_df <- as.data.frame(sub.dat[['RNA']]@data)
expr_df$symbol <- rownames(expr_df)
# have to prepare a gtf 
geneinfo <- read.table("/public/workspace/lily/REF/INDEX-hg19/anno/gencode_hg19_pos.txt")
colnames(geneinfo) <- c("symbol","chr","start","end")
geneinfo <- geneinfo[,c("chr","start","end","symbol")]
geneinfo$chr <- gsub("chr","",geneinfo$chr)

# 1. get gene coord information 
coord.df <- convert_to_coord(expr_df,geneinfo,chrs=1:22)

# 2. normalize 
norm.ge <- norm_ge(coord.df,rcpp = T)

# 3. bin gene
bin.ge <- bin_genes(norm.ge,nb_cores = 4, rcpp = T)

# 4. z score 
zsc.ge <- zscore(bin.ge, normals = Bcell)

# 5. smooth 
smo.ge <- smooth_movingw(zsc.ge, wsize = 3, nb_cores = 4, FUN = stats::median,rcpp = T)

# 6. run PCA 
pca.res <- run_pca(smo.ge, core_cells = NULL, out_pcs = 100)

# 7. find community 
community <- find_communities(pca.res, nb_pcs = 10, k = 100, gamma = 1, nreps = 1,nb_cores = 4)


#============================================================================================
# ready to call CNV 
# 1. make meta cells


# res <- auto_cna_call(ge_df=norm.ge, comm_df=community$comm, nb_metacells = 10)








metainfo <- make_metacells(norm.ge, community$comm, nb_metacells = 10, metacell_size = 3,baseline_cells = Bcell, nb_cores = 4, max_baseline_comm = 3)
meta.ge <- metainfo$ge

# 2. normalize 
norm.meta.ge <- norm_ge(meta.ge,rcpp = T)

# 3. bin gene
bin.meta.ge <- bin_genes(norm.meta.ge,nb_cores = 4, rcpp = T)

# 4. z score 
zsc.meta.ge <- zscore(bin.meta.ge)

# 5. smooth 
smo.meta.ge <- smooth_movingw(zsc.meta.ge, wsize = 3, nb_cores = 4, FUN = stats::median,rcpp = T)

# 6. call cna 
res <- call_cna(smo.meta.ge)





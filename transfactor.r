
#=======================================================================
# 2020-11-3 
#=======================================================================
# analysis baout EMT during tumor metastasis 
#=======================================================================
# maybe not only use rabbit, but also use scenic to find some TFs
#


# Rabbit result 
#=======================================================================
load('/public/workspace/lily/Lung2Brain/inte7/res_Data/TCGA_LUAD_mod.RData') # luad_mod
rabit_luad <- read.table('~/metastasis/data/TCGA_RABIT/LUAD/RABIT_LUAD.HiSeq.V2',sep='\t',header=T) #80 TFs 
rownames(rabit_luad) <- gsub("-",".",rabit_luad[,1])
rabit_luad$X <- NULL		#prepare

# get co-sample 
co_luad <- merge(luad_mod,rabit_luad,by='row.names')
rownames(co_luad) <- co_luad$Row.names
co_luad$Row.names <- NULL
tf_luad <- data.frame(t(apply(co_luad[,-c(1:3)],2,function(x,y){tmp=cor.test(x,y);c(tmp$estimate,tmp$p.value,length(which(x!=0))/nrow(co_luad))},y=co_luad$BMS_test_norm)))
# col 1:3 is mod-result 
colnames(tf_luad) <- c('cor','pvalue','percent')
tf_luad$log2Pvalue <- -log2(tf_luad$pvalue) # -log2(pvalue)
# set color
tf_luad$type <- 'up_regulate'
tf_luad$type[which(tf_luad$cor<0)] <- 'down_regulate'
# filter some TFs
tf_luad.f <- tf_luad[which(tf_luad$pvalue<0.05),]

#====================================================================
# get TFs to use in TRRUST 
up.tf <- rownames(tf_luad.f[which(tf_luad.f$cor>0),])
dn.tf <- rownames(tf_luad.f[which(tf_luad.f$cor<0),])
# write.table(up.tf,file="./tmp_up_tf.txt",quote=F,row.names=F,col.names=F)
# write.table(dn.tf,file="./tmp_dn_tf.txt",quote=F,row.names=F,col.names=F)


trrust.dat <- read.table('/public/workspace/lily/metastasis/data/TRRUST/trrust_rawdata.human.tsv',sep="\t") 
colnames(trrust.dat) <- c("source","traget","effect","id")
# head(trrust.dat[which(trrust.dat$source%in%up.tf),])  # check TFs 
head(trrust.dat[which(trrust.dat$source%in%up.tf&&trrust.dat$effect=="Activation"),])
up.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%up.tf),"traget"]))
dn.gene <-  as.vector(unique(trrust.dat[which(trrust.dat$source%in%dn.tf),"traget"]))

write.table(up.gene,file="./tmp_up_gene.txt",quote=F,row.names=F,col.names=F)
write.table(dn.gene,file="./tmp_dn_gene.txt",quote=F,row.names=F,col.names=F)


#===================================================================
# KEGG analysis 
#
#===================================================================
# 
library()
library(clusterProfiler)
library(org.Hs.eg.db)
geneinfo <- bitr(up.gene, fromType="SYMBOL", toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")	   
ekk <- enrichKEGG(gene= geneinfo$ENTREZID,organism  = 'hsa', qvalueCutoff = 0.05)
a1 <- ekk@result
a1 <- a1[which(a1$pvalue<0.01),]


























































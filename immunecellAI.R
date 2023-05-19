

#=======================================================================================
# caculate TCGA LUAD samples and BMS score 
#
#=======================================================================================
source("./software/ssgseaMOD.r")


#============================ 
# load data
dat <- read.delim("./data/TCGA_LUAD_HiSeqV2",sep="\t",header = T)
head(dat)
rownames(dat) <- dat$sample
dat$sample <- NULL
dim(dat)


mod <- mod.analyze2(as.matrix(dat),c("BMS_test"),"./MOD/",permN=0)
mod <- as.data.frame(mod)
mod$type <- "Median"
mod$type[which(mod$BMS_test_norm>=quantile(mod$BMS_test_norm,0.67))] <- "BMS_high"
mod$type[which(mod$BMS_test_norm<quantile(mod$BMS_test_norm,0.33))] <- "BMS_low"
type <- mod$type
names(type) <- rownames(mod)
dat.f <- rbind(type,dat)
rownames(dat.f)[1] <- "group"
dat.ff <- dat.f[,-which(dat.f[1,]=="Median")]
#=============================
write.table(dat.ff,file="./data/LUAD_BMS_group.txt",sep="\t",quote=F)








#=================================================================================================================================
# immune cell AI result plot 
#=================================================================================================================================
result <- read.delim("./data/ImmuCellAI_icb_result.txt",sep="\t")
rownames(result) <- result$X
res.f <- merge(result,mod,by="row.names")



#=================================================================================================================================
# plot result 
#=================================================================================================================================
boxplot(CD8_naive~type,data=res.f) # CD8_naive
boxplot(Cytotoxic~type,data=res.f) # Cytotoxic 1
boxplot(Exhausted~type,data=res.f) # Exhausted 1
boxplot(Tr1~type,data=res.f) # Tr1
boxplot(nTreg~type,data=res.f) # nTreg
boxplot(iTreg~type,data=res.f) # 
boxplot(Th1~type,data=res.f) # Th1
boxplot(Th2~type,data=res.f) # Th2
boxplot(iTreg~type,data=res.f) # iTreg
boxplot(Macrophage~type,data = res.f) # Macrophage
boxplot(NK~type,data = res.f) # NK
boxplot(CD4_T~type,data=res.f) # CD4_T
boxplot(CD8_T~type,data=res.f) # CD8_T
boxplot(InfiltrationScore~type,data = res.f)









#===================================================================================================================================
# use another GSE data to analysis 
# GSE68465
#===================================================================================================================================
gpl96 <- read.delim("./data/GPL96-57554.txt",sep="\t",comment.char = "#")
gsedat <- read.delim("./data/GSE68465_series_matrix.txt",sep="\t",comment.char = "!",header=T)
gsedat.f <- gsedat[,-(which(colnames(gsedat)=="GSM1672724"):ncol(gsedat))]
tmp <- gpl96[,c("ID","Gene.Symbol")]
colnames(tmp)[1]<- "ID_REF"
gsedat.ff <- merge(tmp,gsedat.f,by="ID_REF")
gsedat.ff$ID_REF <- NULL
fdat.f <- aggregate(.~Gene.Symbol,data=gsedat.ff,FUN=median)

#===================================================================================================================================
source("./software/ssgseaMOD.r")
gsedat <- readRDS("./data/GSE68465_LUAD.RDS")
mod <- mod.analyze2(as.matrix(gsedat),c("BMS_test"),"./MOD/",permN=0)
mod <- as.data.frame(mod)

mod$type <- "Median"
mod$type[which(mod$BMS_test_norm>=quantile(mod$BMS_test_norm,0.67))] <- "BMS_high"
mod$type[which(mod$BMS_test_norm<quantile(mod$BMS_test_norm,0.33))] <- "BMS_low"
type <- mod$type
names(type) <- rownames(mod)
dat.f <- rbind(type,gsedat)
rownames(dat.f)[1] <- "group"
dat.ff <- dat.f[,-which(dat.f[1,]=="Median")]


#=============================
write.table(dat.ff,file="./data/GSE68465_LUAD_BMS_group.txt",sep="\t",quote=F)







#====================================================================================================================================
# GSE135222
#
#====================================================================================================================================
dat <- read.delim("./data/GSE135222_GEO_RNA-seq_omicslab_exp.tsv",sep="\t",header=T)
tmp <- read.delim("./hg38_genelength.txt",sep="\t",header=T)
sapply(strsplit(as.vector(dat$gene_id),"\\."), function(x){x[1]}) -> dat$Gene.stable.ID
dat.f <- merge(tmp[,c(1,4)],dat,by="Gene.stable.ID")
dat.f$Gene.stable.ID <- NULL
dat.f$gene_id <- NULL

#========================================================================================
# select LUAD samples 
#=======================
res.f <- aggregate(.~Gene.name,data = dat.f,FUN=median)
rownames(res.f) <- res.f$Gene.name
res.ff <- res.f[,which(colnames(res.f)%in%paste0("NSCLC",c(378,1203,1401,825,1528,1017,1358,1510,1164,1145,1619,1809,1873)))]


#==================== 
#calculate ssgsea
#=====================
source("./software/ssgseaMOD.r")
mod <- mod.analyze2(as.matrix(res.ff),c("BMS_test"),"./MOD/",permN=0)
mod <- as.data.frame(mod)
mod$type <- "Median"
mod$type[which(mod$BMS_test_norm>=quantile(mod$BMS_test_norm,0.67))] <- "BMS_high"
mod$type[which(mod$BMS_test_norm<quantile(mod$BMS_test_norm,0.33))] <- "BMS_low"
type <- mod$type
names(type) <- rownames(mod)
dat.f <- rbind(type,res.ff)
rownames(dat.f)[1] <- "group"
dat.ff <- dat.f[,-which(dat.f[1,]=="Median")]
write.table(dat.ff,file="./data/GSE135222_LUAD_BMS_group.txt",sep="\t",quote=F)


dat.ff <- dat.f
write.table(dat.ff,file="./data/GSE135222_LUAD_BMS_group.txt",sep="\t",quote=F)




#===================================================================================================================
# plot for immune cell AI result 
#
#===================================================================================================================
#======= TCGA LUAD samples
dat <- read.delim("./data/immuneCellAI/TCGA_LUAD_ImmuCellAI_icb_result.txt",sep="\t",header=T)
colnames(dat)[1] <- "sample"
tmp <- read.delim("./data/LUAD_BMS_group.txt",sep = "\t",header = T)
tmp[1:5,1:5]
typeinfo <- data.frame(sample=colnames(tmp),type=as.vector(unlist(tmp[1,])))
typeinfo$type <- as.vector(typeinfo$type)
typeinfo <- typeinfo[-1,]

res <- merge(dat,typeinfo,by="sample")
res.f <- res
pdf("./data/TCGA_LUAD_immuneCellAI_res.pdf")
boxplot(CD8_naive~type,data=res.f) # CD8_naive
legend("topright",legend = paste0(round(wilcox.test(CD8_naive~type,data=res.f)$p.value,digits=3)))
boxplot(Cytotoxic~type,data=res.f) # Cytotoxic 1
legend("topright",legend = paste0(round(wilcox.test(Cytotoxic~type,data=res.f)$p.value,digits=3)))
boxplot(Exhausted~type,data=res.f) # Exhausted 1
legend("topright",legend = paste0(round(wilcox.test(Exhausted~type,data=res.f)$p.value,digits=3)))
boxplot(Tr1~type,data=res.f) # Tr1
legend("topright",legend = paste0(round(wilcox.test(Tr1~type,data=res.f)$p.value,digits=3)))
boxplot(nTreg~type,data=res.f) # nTreg
legend("topright",legend = paste0(round(wilcox.test(nTreg~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # 
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Th1~type,data=res.f) # Th1
legend("topright",legend = paste0(round(wilcox.test(Th1~type,data=res.f)$p.value,digits=3)))
boxplot(Th2~type,data=res.f) # Th2
legend("topright",legend = paste0(round(wilcox.test(Th2~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # iTreg
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Macrophage~type,data = res.f) # Macrophage
legend("topright",legend = paste0(round(wilcox.test(Macrophage~type,data=res.f)$p.value,digits=3)))
boxplot(NK~type,data = res.f) # NK
legend("topright",legend = paste0(round(wilcox.test(NK~type,data=res.f)$p.value,digits=3)))
boxplot(CD4_T~type,data=res.f) # CD4_T
legend("topright",legend = paste0(round(wilcox.test(CD4_T~type,data=res.f)$p.value,digits=3)))
boxplot(CD8_T~type,data=res.f) # CD8_T
legend("topright",legend = paste0(round(wilcox.test(CD8_T~type,data=res.f)$p.value,digits=3)))
# boxplot(InfiltrationScore~type,data = res.f)
dev.off()












#======================================================================================================
# GSE68546 
#
#=====================================
dat <- read.delim("./data/immuneCellAI/GSE68546_ImmuCellAI_icb_result.txt",sep="\t",header=T)
colnames(dat)[1] <- "sample"
tmp <- read.delim("./data/GSE68465_LUAD_BMS_group.txt",sep = "\t",header = T)
tmp[1:5,1:5]

typeinfo <- data.frame(sample=colnames(tmp),type=as.vector(unlist(tmp[1,])))
typeinfo$type <- as.vector(typeinfo$type)
typeinfo <- typeinfo[-1,]

res <- merge(dat,typeinfo,by="sample")
res.f <- res
pdf("./data/GSE68546_LUAD_immuneCellAI_res.pdf")
boxplot(CD8_naive~type,data=res.f) # CD8_naive
legend("topright",legend = paste0(round(wilcox.test(CD8_naive~type,data=res.f)$p.value,digits=3)))
boxplot(Cytotoxic~type,data=res.f) # Cytotoxic 1
legend("topright",legend = paste0(round(wilcox.test(Cytotoxic~type,data=res.f)$p.value,digits=3)))
boxplot(Exhausted~type,data=res.f) # Exhausted 1
legend("topright",legend = paste0(round(wilcox.test(Exhausted~type,data=res.f)$p.value,digits=3)))
boxplot(Tr1~type,data=res.f) # Tr1
legend("topright",legend = paste0(round(wilcox.test(Tr1~type,data=res.f)$p.value,digits=3)))
boxplot(nTreg~type,data=res.f) # nTreg
legend("topright",legend = paste0(round(wilcox.test(nTreg~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # 
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Th1~type,data=res.f) # Th1
legend("topright",legend = paste0(round(wilcox.test(Th1~type,data=res.f)$p.value,digits=3)))
boxplot(Th2~type,data=res.f) # Th2
legend("topright",legend = paste0(round(wilcox.test(Th2~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # iTreg
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Macrophage~type,data = res.f) # Macrophage
legend("topright",legend = paste0(round(wilcox.test(Macrophage~type,data=res.f)$p.value,digits=3)))
boxplot(NK~type,data = res.f) # NK
legend("topright",legend = paste0(round(wilcox.test(NK~type,data=res.f)$p.value,digits=3)))
boxplot(CD4_T~type,data=res.f) # CD4_T
legend("topright",legend = paste0(round(wilcox.test(CD4_T~type,data=res.f)$p.value,digits=3)))
boxplot(CD8_T~type,data=res.f) # CD8_T
legend("topright",legend = paste0(round(wilcox.test(CD8_T~type,data=res.f)$p.value,digits=3)))
# boxplot(InfiltrationScore~type,data = res.f)
dev.off()







#================================================================================================================
# GSE68546 
#=====================================
#=====================================
dat <- read.delim("./data/immuneCellAI/GSE135222_ImmuCellAI_icb_result.txt",sep="\t",header=T)
colnames(dat)[1] <- "sample"
tmp <- read.delim("./data/GSE135222_LUAD_BMS_group.txt",sep = "\t",header = T)
tmp[1:5,1:5]
typeinfo <- data.frame(sample=colnames(tmp),type=as.vector(unlist(tmp[1,])))
typeinfo$type <- as.vector(typeinfo$type)
typeinfo <- typeinfo[-1,]

res <- merge(dat,typeinfo,by="sample")
res.f <- res
pdf("./data/GSE135222_LUAD_immuneCellAI_res.pdf")
boxplot(CD8_naive~type,data=res.f) # CD8_naive
legend("topright",legend = paste0(round(wilcox.test(CD8_naive~type,data=res.f)$p.value,digits=3)))
boxplot(Cytotoxic~type,data=res.f) # Cytotoxic 1
legend("topright",legend = paste0(round(wilcox.test(Cytotoxic~type,data=res.f)$p.value,digits=3)))
boxplot(Exhausted~type,data=res.f) # Exhausted 1
legend("topright",legend = paste0(round(wilcox.test(Exhausted~type,data=res.f)$p.value,digits=3)))
boxplot(Tr1~type,data=res.f) # Tr1
legend("topright",legend = paste0(round(wilcox.test(Tr1~type,data=res.f)$p.value,digits=3)))
boxplot(nTreg~type,data=res.f) # nTreg
legend("topright",legend = paste0(round(wilcox.test(nTreg~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # 
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Th1~type,data=res.f) # Th1
legend("topright",legend = paste0(round(wilcox.test(Th1~type,data=res.f)$p.value,digits=3)))
boxplot(Th2~type,data=res.f) # Th2
legend("topright",legend = paste0(round(wilcox.test(Th2~type,data=res.f)$p.value,digits=3)))
boxplot(iTreg~type,data=res.f) # iTreg
legend("topright",legend = paste0(round(wilcox.test(iTreg~type,data=res.f)$p.value,digits=3)))
boxplot(Macrophage~type,data = res.f) # Macrophage
legend("topright",legend = paste0(round(wilcox.test(Macrophage~type,data=res.f)$p.value,digits=3)))
boxplot(NK~type,data = res.f) # NK
legend("topright",legend = paste0(round(wilcox.test(NK~type,data=res.f)$p.value,digits=3)))
boxplot(CD4_T~type,data=res.f) # CD4_T
legend("topright",legend = paste0(round(wilcox.test(CD4_T~type,data=res.f)$p.value,digits=3)))
boxplot(CD8_T~type,data=res.f) # CD8_T
legend("topright",legend = paste0(round(wilcox.test(CD8_T~type,data=res.f)$p.value,digits=3)))
# boxplot(InfiltrationScore~type,data = res.f)
dev.off()





































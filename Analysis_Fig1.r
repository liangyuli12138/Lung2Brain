
# Analysis for Fig1 
# this program is used to analysis Figure1 
# 1. landscape - tsne # T cells include NK cells 
# 2. make gene plot 
# 3. cnv correlation with infercnv and TCGA 
# 4. cell type gene marker 

#==============================================================================================================================================
#==============================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
cols <-c('#377EB8','#eeaf00','#f784b6','#984EA3','#A65628','#aea400','#8CA77B')
cols2 <- c('#377EB8','#eeaf00','#f784b6','#984EA3','#A65628','#aea400','#8CA77B','#E41A1C')

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Landscape_tsne.pdf",useDingbats=F,width=8,height=8)
DimPlot(dat,group.by="celltype",reduction="tsne",pt.size=0.35,cols=cols,raster=F)
DimPlot(dat,group.by="celltype.refine",reduction="tsne",pt.size=0.35,cols=cols2,raster=F)
DimPlot(dat,group.by="type_group",reduction="tsne",pt.size=0.35,raster=F,cols=c(rgb(0,160,233,max=255),rgb(243,152,0,max=255),rgb(143,195,31,max=255))) # blue orange green
dev.off()


#==============================================================================================================================================
# 2023-03-04
dat$newgroup <- as.character(dat$type_group)
dat$newgroup[which(dat$orig.ident=="Pair_BM")] <- "Pair_BM"
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Landscape_tsne_split.pdf",useDingbats=F,width=24,height=8)
DimPlot(dat,group.by="celltype.refine",reduction="tsne",split.by="newgroup",pt.size=0.35,cols=cols2,raster=F)
dev.off()



# Bcell Endothelial  Epithelial  Fibroblast     Myeloid      Oligo.
#        6565        1238       20524       18882       38564         767
#       Tcell
#       22949


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Landscape_tsne_type_group.pdf",useDingbats=F,width=8,height=8)
DimPlot(dat,group.by="type_group",reduction="tsne",cols=c("#f7931e","#006837","#6cbc35"))
dev.off()

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Supplementary_Fig1_sample_custers.pdf",useDingbats=F,width=8,height=8)
DimPlot(dat,reduction="tsne",group.by="seurat_clusters")
DimPlot(dat,reduction="tsne",group.by="orig.ident",cols=c('#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1',
'#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#54B0E4','#222F75','#1B9E77'))
dev.off()




#==============================================================================================================================================
#==============================================================================================================================================
# marker gene plot 
# 2022-1-10
library(ggplot2)
marker <- c("CD3D","CD3E","CD2","PTPRC",
			"IGHM","CD79A","MS4A1",
			"CD68","FCGR1A","LYZ",
			"CLDN5","CDH5","VWF",
            "COL1A2","THY1","DCN",
            "MAG","MOG","CNDP1",
			"KRT19","EPCAM","CDH1"
			)
dat$celltype <- factor(dat$celltype,levels=c("Tcell","Bcell","Myeloid","Endothelial","Fibroblast","Oligo.","Epithelial"))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Maker_dotplot.pdf",useDingbats=F)
DotPlot(dat,features=marker,group.by="celltype") + coord_flip() + 
    scale_color_gradientn(colors=rev(c("#FFD92F","#FEE391",RColorBrewer::brewer.pal(11, "Spectral")[7:11]))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



#==============================================================================================================================================
#==============================================================================================================================================
# plot CNV.score between Tumor cells and Epithelial like 
# 2022-3-28
library(Seurat)
library(infercnv)
samplelist<- gsub("\\.RDS$","",grep("*.RDS$",dir("/public/workspace/lily/Lung2Brain/Version6/Prepare_Data/"),value=T))
samplelist <- samplelist[-which(samplelist=="PLCBM2")]
tmp.list1 <- list()
tmp.list2 <- list()
for(i in 1:length(samplelist)){
    respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i],"/")
    dat <- readRDS(paste0(respath,"run.final.infercnv_obj"))
    # calculate all cells CNV.score
    tmp.cnv <- apply(dat@expr.data,2,function(x){sum((x-1)^2)})

    # get reference CNV.score
    # cutoff use 75% of reference CNV score 
    ref.cnv <- quantile(tmp.cnv[dat@reference_grouped_cell_indices$Tcell],0.75)

    # get obs cells 
    obs <- tmp.cnv[dat@observation_grouped_cell_indices$Epithelial]
    obs.res <- obs[which(obs>ref.cnv)]
    print(samplelist[i])
    print(length(dat@observation_grouped_cell_indices$Epithelial))
    print(length(obs.res))
    tmp.list1[[i]] <- obs.res
    names(tmp.list1)[i] <- samplelist[i]
    tmp.list2[[i]] <- obs[-which(obs>ref.cnv)]
    names(tmp.list2)[i] <- samplelist[i]
   
}

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Epithelial_CNV_score.pdf",useDingbats=F)
boxplot(unlist(tmp.list1),unlist(tmp.list2),names=c("Tumor like","Epithelial like"),outline=F,main="CNV.score(sum of squares)")
# boxplot(sqrt(unlist(tmp.list1)),sqrt(unlist(tmp.list2)),names=c("Tumor like","Epithelial like"),outline=F,main="CNV.score(sum of squares)")
dev.off()



# 2022-5-18 
# another plot show Tumor in differnet group with Epithelial
library(Seurat)
library(ggplot2)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
cnv.score <- unlist(c(tmp.list1,tmp.list2))
names(cnv.score) <- unique(gsub("^.*\\.","",names(unlist(cnv.score))))
cnv.score <- cnv.score[colnames(dat)]
dat$cnv.score <- cnv.score
dat$tmp.type <- paste0(dat$type_group,"_",dat$celltype.refine)
dat$tmp.type[grep("Epithelial",dat$tmp.type)] <- "Epithelial"


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Epithelial_CNV_score_group.pdf",useDingbats=F)
boxplot(cnv.score~tmp.type,data=dat@meta.data,outline=F,main="CNV.score(sum of squares)")
# boxplot(sqrt(unlist(tmp.list1)),sqrt(unlist(tmp.list2)),names=c("Tumor like","Epithelial like"),outline=F,main="CNV.score(sum of squares)")
dev.off()




# and add cell num plot 
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")
dat$tmp.type <- paste0(dat$type_group,"_",dat$celltype.refine)
dat$tmp.type[grep("Epithelial",dat$tmp.type)] <- "Epithelial"

library(plotrix)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Epithelial_cellnum.pdf",useDingbats=F)
barplot(table(dat$tmp.type),names=names(table(dat$tmp.type)),ylim=c(0,5000))
dev.off()

























#==============================================================================================================================================
#==============================================================================================================================================
# now run CNV for 3 samples
library(infercnv)
library(Seurat)

scCNV2chr <- function(sampleid){
    dat <- readRDS(paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",sampleid,"/run.final.infercnv_obj"))
    expr <- as.data.frame(dat@expr.data[,-dat@reference_grouped_cell_indices$Tcell])
    expr$gene_order <- dat@gene_order$chr
    res <- aggregate(.~gene_order,data=expr,FUN=mean) #get every chr for each cell 
    rownames(res) <- res$gene_order
    res$gene_order <- NULL
    res.f <- apply(res,1,function(x){mean(x)})
    return(res.f)
}

segCNV2chr <- function(dat){
    # make sure your colnames
    dat <- data.frame(dat)
    dat$len <- dat$loc.end -dat$loc.start
    dat$seg.mean <- as.numeric(as.vector(dat$seg.mean))
    dat$score <- dat$len*dat$seg.mean
    # need one more step to calculate values 
    res.score <- aggregate(score~chrom,data=dat,FUN=sum)
    res.len <- aggregate(len~chrom,data=dat,FUN=sum)
    res.f <- merge(res.score,res.len,by="chrom")
    res.f$value <- res.f$score/res.f$len
    return(res.f)

}

#===============================================================================================================================================
# read CNVkit information 

tmp <- read.table("/public/workspace/lily/Lung2Brain/Version6/CNV/CNVkit/A20190305_T.seg",header=T,sep="\t")
scCNV.A05 <- scCNV2chr("A20190305")
scCNV.A05 <- scCNV.A05[paste0("chr",1:22)]
seg.A05 <- segCNV2chr(tmp)
rownames(seg.A05) <- seg.A05$chr
seg.A05 <- seg.A05[paste0("chr",1:22),]
cor.test(seg.A05$value,(scCNV.A05-1),method="spearman")

tmp <- read.table("/public/workspace/lily/Lung2Brain/Version6/CNV/CNVkit/A20190312_T.seg",header=T,sep="\t")
scCNV.A12 <- scCNV2chr("A20190312")
scCNV.A12 <- scCNV.A12[paste0("chr",1:22)]
seg.A12 <- segCNV2chr(tmp)
rownames(seg.A12) <- seg.A12$chr
seg.A12 <- seg.A12[paste0("chr",1:22),]
cor.test(seg.A12$value,(scCNV.A12-1),method="spearman")

tmp <- read.table("/public/workspace/lily/Lung2Brain/Version6/CNV/CNVkit/T_Bsc_T.seg",header=T,sep="\t")
scCNV.TBsc <- scCNV2chr("T_Bsc1")
scCNV.TBsc <- scCNV.TBsc[paste0("chr",1:22)]
seg.TBsc <- segCNV2chr(tmp)
rownames(seg.TBsc) <- seg.TBsc$chr
seg.TBsc <- seg.TBsc[paste0("chr",1:22),]
cor.test(seg.TBsc$value,(scCNV.TBsc-1),method="spearman")


# plot result

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/WES_CNV_scCNV_cor.pdf",useDingbats=F)
plot(seg.A12$value,(scCNV.A12-1),main="A12 low tumor cells")
abline(lm((scCNV.A12-1)~seg.A12$value),col="red")
legend("topright",legend=paste0("rho=",cor.test(seg.A12$value,(scCNV.A12-1),method="spearman")$estimate,
    " P =",round(cor.test(seg.A12$value,(scCNV.A12-1),method="spearman")$p.value,2))
)

plot(seg.TBsc$value,(scCNV.TBsc-1),main="T_Bsc")
abline(lm((scCNV.TBsc-1)~seg.TBsc$value),col="red")
legend("topright",legend=paste0("rho=",cor.test(seg.TBsc$value,(scCNV.TBsc-1),method="spearman")$estimate,
    " P =",round(cor.test(seg.TBsc$value,(scCNV.TBsc-1),method="spearman")$p.value,2))
)

plot(seg.A05$value,(scCNV.A05-1),main="A05")
abline(lm((scCNV.A05-1)~seg.A05$value),col="red")
legend("topright",legend=paste0("rho=",cor.test(seg.A05$value,(scCNV.A05-1),method="spearman")$estimate,
    " P =",round(cor.test(seg.A05$value,(scCNV.A05-1),method="spearman")$p.value,2))
)

dev.off()









#==============================================================================================================================================
#==============================================================================================================================================
# now plot CNV for 16 samples
# bytlib load JAGS-4.3.0
library(Seurat)
library(infercnv)


CalCulateCNV=function(infercnv_obj,hpdat,SampleName,chrinfo,mean_hpdat){
 options(stringsAsFactors=F)
 library(infercnv,lib.loc='/public/workspace/lily/R/x86_64-pc-linux-gnu-library/3.6.0/')
 library(TxDb.Hsapiens.UCSC.hg38.knownGene)
 library(ggplot2)

 DisPlayGene=infercnv_obj@gene_order
 ChromosomeList=paste0('chr',1:22)
 chrinfo=subset(chrinfo,chrom%in%ChromosomeList)
 DisPlayGene=subset(DisPlayGene,chr%in%ChromosomeList)
 chr.gr = GRanges(seqnames=chrinfo$chrom,ranges=IRanges(start=chrinfo$chromStart,end=chrinfo$chromEnd))
 cnv.gr = GRanges(seqnames=DisPlayGene$chr,ranges=IRanges(start=DisPlayGene$start,end=DisPlayGene$stop),mcols=as.numeric(signif(mean_hpdat,3)))
 ind = findOverlaps(cnv.gr,chr.gr)
 pos=chrinfo[subjectHits(ind),c('chrom','name')]
 cnv.gr=cnv.gr[queryHits(ind),]
 mcols(cnv.gr)$pos=pos$name
 mcols(cnv.gr)$chr=pos$chrom
 cnv.gr=as.data.frame(cnv.gr)
 cnv.gr=unique(cnv.gr)
 seg=as.data.frame(cnv.gr%>%group_by(pos,chr)%>%summarise(seg=mean(mcols)))
 seg$chr=factor(seg$chr,levels=rev(paste0('chr',1:22)))
 seg$pos=factor(seg$pos,levels=c('q','p'))
 seg$Sample=SampleName
 return (seg)
}





library(data.table)
library(ggplot2)
library(dplyr)
chrinfo <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
chrinfo$name=gsub('[0-9].*','',chrinfo$name)
ChromosomeList=paste0('chr',1:22)

objfile=Sys.glob('/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/*/run.final.infercnv_obj')
obvfile=Sys.glob('/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/*/infercnv.observations.txt')

#=========================================================================================================
SEG=data.frame()
for(i in 1:length(objfile)){
    infercnv_obj=readRDS(objfile[i]);
    hpdat=read.table(obvfile[i],header=T,row.names=1,sep=' ',stringsAsFactors=F)
    # 2022-6-27 modify just use tumor cells
    tumor <- readRDS("~/Lung2Brain/Version6/Data/Tumor.RDS")
    hpdat <- hpdat[,which(colnames(hpdat)%in%colnames(tumor))]
    cat(length(which(colnames(hpdat)%in%colnames(tumor))))
    cat("\n")
    hpdat=signif(hpdat,3)
    # mean_hpdat=apply(hpdat,1,mean)
    # highvalue <- mean(mean(as.numeric(hpdat["ACTB",])),mean(as.numeric(hpdat["GAPDH",]))) + 0.1
    # lowvalue <- highvalue - 0.2
    mean_hpdat=apply(hpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})

    mean_hpdat=mean_hpdat/ncol(hpdat)
    SampleName=gsub("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/|/infercnv.observations.txt","",obvfile[i])
    SEG=rbind(SEG,CalCulateCNV(infercnv_obj, hpdat,SampleName,chrinfo,mean_hpdat))
}
targetinfo=data.frame(chr=rep(rep(ChromosomeList,rep(2,length(ChromosomeList))),length(unique(SEG$Sample))),
pos=rep(c('p','q'),length(ChromosomeList)*length(unique(SEG$Sample))),
Sample=rep(unique(SEG$Sample),rep(2*length(ChromosomeList),length(unique(SEG$Sample)))))
SEG=dplyr::left_join(targetinfo,SEG,by=c('chr','pos','Sample'))
SEG[is.na(SEG$seg),'seg']=0
SEG$chr=factor(SEG$chr,levels=rev(paste0('chr',1:22)))
SEG$pos=factor(SEG$pos,levels=c('q','p'))

SEG$Sample <- factor(SEG$Sample,levels=c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T25","LUNG_T30","LUNG_T34",
    "Pair_Lung","Pair_BM","A20190305","A20190312","D0927","E0927","T_Bsc1"))



# plot function
plotCNV=function(SEG,limit=c(-0.2,0.2)){
 library(ggplot2)
 idx = rev(table(unique(SEG[,c(1,2)])$chr))
 idx[idx>=2]=2
 height=c(0)
 for(i in 1:length(idx)){height=c(height,sum(idx[1:i]))}
 bk=c()
 for(i in 1:length(idx)){bk=c(bk,sum(idx[1:i])-idx[i])}
 SEG = arrange(SEG,Sample,desc(chr),desc(pos))
 MIN=(limit[1]-min(SEG$seg))/(max(SEG$seg)-min(SEG$seg))
 MAX=(limit[2]-min(SEG$seg))/(max(SEG$seg)-min(SEG$seg))
 
 p=ggplot() + 
  geom_bar(data = SEG, aes(x = Sample,y=1,fill = seg,group=chr),stat = "identity", lwd = 0,color='NA')
  p=p+geom_hline(aes(yintercept=height), colour="black",lwd=0.5)
  p=p+coord_flip()+scale_fill_gradientn(colours = c('steelblue',"white","white",'darkred'),
  values=c(0,MIN,MAX,1))
  p=p+scale_y_continuous(name = "chromosome",breaks = bk,label=rev(levels(SEG$chr)),expand=c(0,0))+
  scale_x_discrete(name="",expand=c(0,0))
  p=p+ theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))
  p=p+theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),
   panel.grid.minor.y = element_blank())
  p=p+theme(axis.text.x = element_text(angle=45,size=7,hjust=1))
  p=p+theme(legend.position="bottom")
 return (p)
}

p=plotCNV(SEG)
p

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/CNV_replot_sample.pdf",useDingbats=F,width=8,height=6)
p
dev.off()

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/CNV_replot_sample_tumor.pdf",useDingbats=F,width=8,height=6)
p
dev.off()

pdf("~/tmp.pdf",useDingbats=F)
samplelist <- unique(SEG$Sample)
for(i in 1:length(samplelist)){
   lims <- c(as.numeric(quantile(SEG[which(SEG$Sample==samplelist[i]),"seg"],0.3)),as.numeric(quantile(SEG[which(SEG$Sample==samplelist[i]),"seg"],0.6)))
   p <- plotCNV(SEG[which(SEG$Sample==samplelist[i]),],limit=lims)
   print(p)
}
dev.off()






#==============================================================================================================================================
#==============================================================================================================================================
#=======================================================================================================================================
# calculate signifcant for CNV
# 2022-5-20

tmp <- SEG
tmp$type <- paste0(tmp$chr,tmp$pos)
tmp$group <- "MLUAD"

tmp$group[which(tmp$Sample%in%c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T25","LUNG_T30","LUNG_T34"))] <- "nMLUAD"
tmp$group[which(tmp$Sample%in%c("Pair_BM","A20190305","A20190312","D0927","E0927","T_Bsc1"))] <- "LCBM"
tmp <- tmp[-grep("^MLUAD",tmp$group),] # do not calculate this group only 1 sample

for(i in unique(tmp$type)){   
    tmp.dat <- tmp[which(tmp$type == i),]
    print(i)
    print(wilcox.test(seg~group,data=tmp.dat)$p.value)
}

# calculate significant 
# chr1q,chr5p, chr21q
# plot result 
pdat <- tmp[which(tmp$type%in%c("chr1q","chr5p","chr21q")),]
pdat$group <- factor(pdat$group,levels=c("nMLUAD","LCBM"))
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/CNV_replot_chrpq_signifcant.pdf",useDingbats=F,width=8,height=6)

ggplot(pdat,aes(x=type,y=seg,fill=group))+
	geom_boxplot(width=0.5)+
	geom_point(data = pdat, aes(x=type,y = seg,fill=group),position = position_jitterdodge(),size = 3, shape = 21)+
	scale_fill_manual(values = c("#d4c78c","#5ec6f2","#003468","#4b1702"))+
    # facet_wrap(~type) + 
    theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))

dev.off()

#=======================================================================================================================================
#=======================================================================================================================================
# not O.K.
# use BMS score to do correlation
# 2022-9-23 # use Chr5p as group is OK
# TCGA segment check chr5p
tmp.dat <- read.table("/public/workspace/lily/Lung2Brain/Data/luad_tcga_pub_segments.seg",header=T,sep="\t")

dat <- tmp.dat[which(tmp.dat$chrom==5&tmp.dat$loc.start<46100000&tmp.dat$loc.end<46100000),]
dat$length <- dat$loc.end - dat$loc.start
dat$value <- dat$seg.mean * dat$length
tmp.res <- aggregate(.~ID,data=dat[,c(1,7,8)],FUN=sum)
tmp.res$res <- tmp.res$value / tmp.res$length
tmp.res$ID <- gsub("-",".",tmp.res$ID)
tmp.res$status <- ifelse(tmp.res$res>0,"gain","loss")
# load BMS update score 
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
mod$ID <- rownames(mod)
exp <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
all(colnames(exp)==rownames(mod))
mod$OSMR <- as.numeric(exp["OSMR",])
mod$OSM <- as.numeric(exp["OSM",])

res <- merge(tmp.res,mod[,-c(1,3)],by="ID")

#
res$group <- "none"
res$group[which(res$res<0)] <-"loss"
res$group[which(res$res>0.5)] <-"gain"

wilcox.test(OSMR~group,data=res[which(res$group%in%c("gain","loss")),],FUN=median)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/TCGA_LUAD_CHr5p_OSMR.pdf",useDingbats=F)
boxplot(OSMR~group,data=res)
beeswarm::beeswarm(OSMR~group,data=res[which(res$group%in%c("gain","loss")),],col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

boxplot(OSM~group,data=res)
beeswarm::beeswarm(OSM~group,data=res[which(res$group%in%c("gain","loss")),],col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

boxplot(BMS_V6_norm~group,data=res)
beeswarm::beeswarm(BMS_V6_norm~group,data=res[which(res$group%in%c("gain","loss")),],col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

dev.off()

# res <- readRDS("~/tmp/chr5p.RDS") # from MSGDB
# make some files about chr5p
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# for(i in 1:length(res)){
#    mod.generate(res[[i]],names(res)[i],out=paste0("/public/workspace/lily/MOD_file/",names(res)[i],".mod")) # make a mod file  

# }
# mod.generate(unique(unlist(res)),"Chr5p",out=paste0("/public/workspace/lily/MOD_file/Chr5p.mod"))

# dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),names(res),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)

# score <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
# mod$BMS <- score[,2]

# use mod to calculate correlation



#=======================================================================================================================================
#=======================================================================================================================================
# use this data to calculate score 
# GSE200563
#=======================================================================================================================================
# 2022-5-21
# plot result 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
tmp.res <- readRDS("~/metastasis/data/verify/GSE200563/TCGA_stage_I_II_GSE200563_BM_MLUAD_list.RDS")

mod.gse <- as.data.frame(mod.analyze2(as.matrix(tmp.res$GSE_Data),"Chr5p","/public/workspace/lily/MOD_file/",permN=0))
mod.tcga <- as.data.frame(mod.analyze2(as.matrix(tmp.res$TCGA_Data),"Chr5p","/public/workspace/lily/MOD_file/",permN=0))


plotdat <- data.frame(group=c(rep("TCGA",nrow(mod.tcga)),as.vector(tmp.res$GSE_info$group)),score.chr5p=c(mod.tcga[,2],mod.gse[,2]))

# make a data frame to plot violin 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/TCGA_GSE200563_chr5p.pdf",useDingbats=F)
boxplot(score.chr5p~group,plotdat,outline=F)
dev.off()





























# check chr5p gene with BMS score correlation
#============================================================================================================================
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
score <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")
tmp.gene <- unique(unlist(readRDS("~/tmp/chr5p.RDS")))
gene <- tmp.gene[which(tmp.gene%in%rownames(dat))] # filter some gene

tmp.res <- sort(sapply(gene,function(x){
    cor.test(as.numeric(as.vector(dat[x,])), score[,2],method="spearman" )$estimate  
}))

# OSMR.rho   TRIP13.rho
# 0.526818888  0.546319915











#==============================================================================================================================================
#==============================================================================================================================================
# calculate cell stemness and cell cycle
library(Seurat)
library(CytoTRACE)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")

# cell stemness cytoTrace
results <- CytoTRACE(as.matrix(dat[["RNA"]]@data),ncores = 4,subsamplesize = 1000)
dat$cytotrace <- results$CytoTRACE

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Tumor_cytotrace.pdf",useDingbats=F)
boxplot(cytotrace~type_group,data=dat@meta.data,outline=F)
dev.off()



# cell cycle state 
dat <- CellCycleScoring(dat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)



















#==============================================================================================================================================
#==============================================================================================================================================
# cell type percentage plot 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

# and show cell type percentage 
tmp <- table(dat$type_group,dat$celltype)
tmp.f <- tmp[,-grep("Epithelial",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("nMLUAD","MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Endothelial","Fibroblast","Myeloid","Oligo.","Tcell"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Ntumor_celltype.pdf",useDingbats=F)
cols <- c('#377EB8','#910241','#984EA3','#F29403','#aea400','#B2DF8A')
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()


# 2023-02-25
# plot by bar plot 
tmp <- table(dat$orig.ident,dat$celltype)
tmp.f <- tmp[,-grep("Epithelial",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

library(reshape)
library(ggplot2)
library(ggalluvial)
library(ggbreak)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Group <- "nMLUAD"
tmp.dat$Group[which(tmp.dat$Samples%in%c("Pair_Lung"))] <- "MLUAD"
tmp.dat$Group[which(tmp.dat$Samples%in%c("T_Bsc1","Pair_BM","D0927","E0927","A20190312","A20190305"))] <- "LCBM"
tmp.dat$Group <- factor(tmp.dat$Group,levels=c("nMLUAD","MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Endothelial","Fibroblast","Myeloid","Oligo.","Tcell"))
tmp.dat$Group2 <- as.character(tmp.dat$Group)
tmp.dat$Group2[which(tmp.dat$Sample%in%"Pair_Lung")] <- "Pair"
tmp.dat$Group2[which(tmp.dat$Sample%in%"Pair_BM")] <- "Pair"
tmp.dat$Group2 <- factor(tmp.dat$Group2,levels=c("nMLUAD","Pair","LCBM"))



pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Ntumor_celltype_barplot.pdf",useDingbats=F)

ggplot(tmp.dat,aes(x=Cell_type,y=value,fill=Group))+
	geom_bar(fun="median",position=position_dodge(0.8),stat="summary",width=0.5,size=1)+
	stat_summary(fun.data = 'mean_se', geom = "errorbar", colour = "black",width = 0.2,position = position_dodge(0.8))+
	geom_point(data = tmp.dat, aes(x=Cell_type,y = value,fill=Group2),position = position_dodge(0.8),size = 3, shape = 21)+
	scale_fill_manual(values = c(rgb(0,160,233,max=255),rgb(243,152,0,max=255),rgb(143,195,31,max=255),"red"))+ theme_bw() 
# +
	# theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+scale_y_break(c(0.5,0.75),#截断位置及范围
    #             space = 0.3,#间距大小
    #             scales = 0.8)#上下显示比例，大于1上面比例大，小于1下面比例大

dev.off()

# calculate significant 
# by fisher

rs.list <- list()
types <- unique(tmp.dat$Cell_type)
for(i in 1:length(types)){
   rs.list[[i]] <-  c(
    wilcox.test(value~Group,data=tmp.dat[which(tmp.dat$Cell_type==as.character(types[i])&tmp.dat$Group%in%c("LCBM","nMLUAD")),])$p.value # MLUAD and nMLUAD
    )
    names(rs.list)[i] <- as.character(types[i])
}

# $Bcell
# [1] 0.06633367

# $Endothelial
# [1] 0.8639361

# $Fibroblast
# [1] 0.02557443

# $Myeloid
# [1] 0.02557443

# $Oligo.
# [1] 0.0004262389

# $Tcell
# [1] 0.0007992008



#  Endothelial  Epithelial  Fibroblast     Myeloid      Oligo.
#        6565        1238       20524       18882       38564         767
#       Tcell
#       22949


# another ways to plot result 
#==============================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

# and show cell type percentage 
tmp <- table(dat$orig.ident,dat$celltype)
tmp.f <- tmp[,-grep("Epithelial|Endothelial|Fibroblast",colnames(tmp))]
res.f <- apply(tmp.f,1,function(x){x/sum(x)}) # 

res.f <- t(apply(res.f,1,function(x){
    c(median(x[c(1:4,14,16)]),median(x[5:13]),median(x[15]))
}))
colnames(res.f) <- c("LCBM","nMLUAD","MLUAD")

library(reshape)
library(ggplot2)
library(ggalluvial)
tmp.dat <- melt(res.f,id="col.names")
colnames(tmp.dat) <- c("Cell_type","Samples","value")
tmp.dat$Samples <- factor(tmp.dat$Samples,levels=c("nMLUAD","MLUAD","LCBM"))
tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Endothelial","Fibroblast","Myeloid","Oligo.","Tcell"))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Ntumor_celltype.pdf",useDingbats=F)
cols <- c('#377EB8','#910241','#984EA3','#F29403','#aea400','#B2DF8A')
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



tmp.dat$Cell_type <- factor(tmp.dat$Cell_type,level=c("Bcell","Myeloid","Oligo.","Tcell"))
cols <- c('#377EB8','#F29403','#aea400','#B2DF8A')
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/Immune_celltype.pdf",useDingbats=F)
ggplot(tmp.dat, aes(x = Samples, y = value, fill = Cell_type,stratum = Cell_type, alluvium = Cell_type)) +
geom_stratum(width=0.45) +  
geom_flow(alpha = 0.5) +
scale_fill_manual(values = cols)+
theme_bw() + theme(panel.grid.major = element_blank())
dev.off()



# calculate significant 
# by fisher
tmp <- table(dat$type_group,dat$celltype)
rs.list <- list()
for(i in 1:ncol(tmp)){
   rs.list[[i]] <-  c(
    fisher.test(matrix(c(tmp[2,i],(rowSums(tmp)[2]-tmp[2,i]),tmp[2,i],(rowSums(tmp)[3]-tmp[3,i])),ncol=2))$p.value, # MLUAD and nMLUAD
    fisher.test(matrix(c(tmp[1,i],(rowSums(tmp)[1]-tmp[1,i]),tmp[3,i],(rowSums(tmp)[3]-tmp[3,i])),ncol=2))$p.value # LCBM and nMLUAD
    )
    names(rs.list)[i] <- colnames(tmp)[i]
}







##################################################################################################################################################
# 2023-02-25
# use estimate to calculate 
# 用estimate算，感觉BUlk数据和取样还有纯度关系很大，都是在BM中基质浸润，免疫浸润更低
# library(estimate)
# tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
# tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# # just filter MLUNG samples
# sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM","MLUNG")),]
# dat <- tmp.dat[,rownames(sampleinfo)]

# result1 <- CIBERSORT('~/software/LM22.txt','~/software/GSE14108_exp.txt', perm = 100, QN = T)
# result2 <- CIBERSORT('~/software/LM22.txt','~/metastasis/data/verify/TCGA_LUAD/LUAD_RNAseq_Exp.txt', perm = 100, QN = T)
# result3 <- CIBERSORT('~/software/LM22.txt','~/software/GSE200563_MLUNG_BM_exp.txt', perm = 100, QN = T)

# library(estimate)
# filterCommonGenes(input.f="~/software/GSE200563_MLUNG_BM_exp.txt", 
#     output.f="/public/workspace/lily/metastasis/data/verify/GSE200563/GSE200563_MLUNG_BM.gct", id="GeneSymbol")
# estimateScore(input.ds="/public/workspace/lily/metastasis/data/verify/GSE200563/GSE200563_MLUNG_BM.gct",
#     output.ds="/public/workspace/lily/metastasis/data/verify/GSE200563/GSE200563_MLUNG_BM_estimate_score.gct", platform="illumina")

# # GSE4108
# library(estimate)
# filterCommonGenes(input.f="~/software/GSE14108_exp.txt", 
#     output.f="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108.gct", id="GeneSymbol")
# estimateScore(input.ds="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108.gct",
#     output.ds="/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_estimate_score.gct", platform="illumina")



# LUAD_est <- read.table("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/LUAD_estimate_score.gct", skip = 2, header = TRUE)
# GSE200563_est <- read.table("/public/workspace/lily/metastasis/data/verify/GSE200563/GSE200563_MLUNG_BM_estimate_score.gct", skip = 2, header = TRUE)
# GSE14108_est <- read.table("/public/workspace/lily/metastasis/data/verify/GSE14108/GSE14108_estimate_score.gct", skip = 2, header = TRUE)


# pdat <- data.frame(Immunescore= c(as.numeric(LUAD_est[2,3:ncol(LUAD_est)]),as.numeric(GSE200563_est[2,3:ncol(GSE200563_est)])),
#     Stromamscore=c(as.numeric(LUAD_est[1,3:ncol(LUAD_est)]),as.numeric(GSE200563_est[1,3:ncol(GSE200563_est)])),
#     Group=("LUAD"),stringsAsFactors=F)
# pdat$Sample <- c(colnames(LUAD_est[2,3:ncol(LUAD_est)]),colnames(GSE200563_est[2,3:ncol(GSE200563_est)]))
# pdat$Group[which(pdat$Sample%in%rownames(sampleinfo)[which(sampleinfo$group=="MLUNG")])] <- "MLUNG"
# pdat$Group[which(pdat$Sample%in%rownames(sampleinfo)[which(sampleinfo$group=="BM")])] <- "BM"


# pdat <- data.frame(Immunescore= c(as.numeric(LUAD_est[2,3:ncol(LUAD_est)]),as.numeric(GSE14108_est[2,3:ncol(GSE14108_est)])),
#     Stromamscore=c(as.numeric(LUAD_est[1,3:ncol(LUAD_est)]),as.numeric(GSE14108_est[1,3:ncol(GSE14108_est)])),
#     Group=("LUAD"),stringsAsFactors=F)
# pdat$Sample <- c(colnames(LUAD_est[2,3:ncol(LUAD_est)]),colnames(GSE14108_est[2,3:ncol(GSE14108_est)]))
# pdat$Group[which(pdat$Sample%in%colnames(GSE14108_est))] <- "BM"






# use cibersort result 
#=================================================================================================================================================
# 2023-02-25
# from Figure T cells Fig 6
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("BM","MLUNG")),]
dat <- tmp.dat[,rownames(sampleinfo)]
# write.table(dat,file="~/software/GSE200563_MLUNG_BM_exp.txt",sep="\t",row.names=T,col.names=T,quote=F)

LUAD <- readRDS("~/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_cibersort.RDS")
GSE200563_cib <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_cibsersort.RDS")
# saveRDS(result1,file="~/metastasis/data/verify/GSE14108/GSE14108_cibsersort.RDS")

tmp.dat <- rbind(LUAD,GSE200563_cib)
pdat <- data.frame(reshape2::melt(tmp.dat),stringsAsFactors=F)
colnames(pdat) <- c("Sample","Type","Value")
pdat$Group <- "LUAD"
pdat$Group[which(pdat$Sample%in%rownames(sampleinfo)[which(sampleinfo$group=="MLUNG")])] <- "MLUNG"
pdat$Group[which(pdat$Sample%in%rownames(sampleinfo)[which(sampleinfo$group=="BM")])] <- "BM"

pdat$Group <- factor(pdat$Group,levels=c("LUAD","MLUNG","BM"))
library(ggplot2)
library(ggpubr)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/GSE200563_LUAD_cibersort_boxplot.pdf",useDingbats=F,width=30,height=9)
ggplot(pdat,aes(x=Type,y=Value,fill=Group))+
	geom_boxplot(width=1,outlier.shape = NA)+
	# geom_point(data = pdat, aes(x=type,y = seg,fill=group),position = position_jitterdodge(),size = 3, shape = 21)+
	scale_fill_manual(values = c(rgb(143,195,31,max=255),rgb(243,152,0,max=255),rgb(0,160,233,max=255)))+
    theme_bw() +
	theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))+
    stat_compare_means(aes(group=Group),method="anova")

dev.off()
































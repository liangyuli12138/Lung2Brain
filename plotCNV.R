#cnv=read.table('C:/Users/wangzy/Desktop/LD7f7033_segments.txt',comment.char='#',sep='\t',stringsAsFactors=F,header=T);
#####################################################
############ plot the result of sequenza ############
############ @cnv segment.txt including chr start end segmean
############ @GenomeVersion hg19/hg37/hg38
############ @ChromosomeList the chromosome you want to show
############ @ChromWidth Whether the chromosomes are of equal width
############ @limit normal ranges
plotGemomic=function(cnv,limit=c(2,2),GenomeVersion='hg19',ChromosomeList=paste0('chr',c(1:22,'X')),
	maxvalue=4,minvalue=0,ChromWidth=TRUE){
	library(dplyr)
	library(reshape2)
	library(ggplot2)
	if('ALL'%in%ChromosomeList){ChromosomeList=paste0('chr',c(1:22,'X'))}
	if(length(grep('chr',ChromosomeList))==0){ChromosomeList=paste0('chr',ChromosomeList)}
#STEP1. 判断基因组版本，获取选定染色体长度，并创建GRanges对象
	if(GenomeVersion%in%c('hg19','hg37')){
		library(TxDb.Hsapiens.UCSC.hg19.knownGene)
		LENGTH=seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)
	}
	if(GenomeVersion%in%c('hg38')){
		library(TxDb.Hsapiens.UCSC.hg38.knownGene)
		LENGTH=seqlengths(TxDb.Hsapiens.UCSC.hg38.knownGene)
	}
	LENGTH=LENGTH[ChromosomeList]
	#cnv$chromosome=paste0('chr',cnv$chromosome)
	cnv=subset(cnv,chromosome%in%ChromosomeList)
	
	chrinfo=data.frame(
		chromosome=names(LENGTH),
		start.pos=1,
		end.pos=as.numeric(LENGTH)
	)

	chr.gr = GRanges(
		seqnames=chrinfo$chromosome,
		ranges=IRanges(start=chrinfo$start.pos,end=chrinfo$end.pos))
#STEP2. 依据CNV,创建Granges对象	
	cnv.gr = GRanges(
		seqnames=cnv$chromosome,
		ranges=IRanges(start=cnv$start.pos,end=cnv$end.pos),
		segmean= cnv$CNt
		)
#STEP3. 获取差集，填补为正常	
	Gap.gr = setdiff(chr.gr,cnv.gr)
	Gap.gr = GRanges(
		seqnames=Gap.gr@seqnames,
		ranges=Gap.gr@ranges,
		segmean= rep(mean(limit),length(Gap.gr))
		)
	cnv.gr=sort(c(Gap.gr,cnv.gr))
	cnv.gr=as.data.frame(cnv.gr)
	cnv.gr$width=cnv.gr$end-cnv.gr$start+1

	#maxvalue=max(limit)+(max(limit)-min(cnv.gr$segmean))
	cnv.gr[cnv.gr$segmean>maxvalue,'segmean']=maxvalue
	cnv.gr$Sample='Sample1'
	cnv.gr$seqnames=factor(cnv.gr$seqnames,levels=rev(ChromosomeList))
	
#STEP4. 画图
	if(ChromWidth){
		if(length(ChromosomeList)==1){
			height=LENGTH[1];bk=height[1]/2
		}else{
			height=sapply(1:(length(LENGTH)-1),function(x){sum(LENGTH[1:x])})
			bk=c(height[1]/2,sapply(2:length(height),function(x){mean(height[(x-1):x])}),(height[length(height)]+(LENGTH[length(LENGTH)]/2)))
		}
		p=ggplot() + 
		geom_bar(data = cnv.gr, aes(x = Sample,y=width,fill = segmean,group=seqnames),stat = "identity", lwd = 0,color='NA')
		p=p+geom_hline(aes(yintercept=height), colour="black",lwd=0.5)
		p=p+coord_flip()+scale_fill_gradient2(low = "steelblue", mid = "white",high = "darkred",midpoint=mean(limit))
		p=p+scale_y_continuous(name = "chromosome",breaks = bk,label=rev(levels(cnv.gr$seqnames)),expand=c(0,0))+
		scale_x_discrete(name="",breaks = NULL,expand=c(0,0))
		p=p+ theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))
		p=p+theme(axis.text.x = element_text(angle=45,size=7,hjust=1))
		p=p+theme(legend.position="bottom")
	}else{	
		LENGTH=data.frame(seqnames=names(LENGTH),ChromWidth=as.numeric(LENGTH))
		cnv.gr=dplyr::inner_join(cnv.gr,LENGTH,by = "seqnames")
		cnv.gr$WidthRatio=cnv.gr$width/cnv.gr$ChromWidth
		cnv.gr$seqnames=factor(cnv.gr$seqnames,levels=rev(ChromosomeList))
		height=1:(nrow(LENGTH)-1)
		bk=seq(0.5,nrow(LENGTH),1)
		p=ggplot() + 
		geom_bar(data = cnv.gr, aes(x = Sample, y=WidthRatio,fill = segmean,group=seqnames),stat = "identity", lwd = 0,color='NA')
		p=p+geom_hline(aes(yintercept=height), colour="black",lwd=0.5)
		p=p+coord_flip()+scale_fill_gradient2(low = "steelblue", mid = "white",high = "darkred",midpoint=mean(limit))
		p=p+scale_y_continuous(name = "chromosome",breaks = bk,label=rev(levels(cnv.gr$seqnames)),expand=c(0,0))+
		scale_x_discrete(name="",breaks = NULL,expand=c(0,0))
		p=p+ theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))
		p=p+theme(axis.text.x = element_text(angle=45,size=7,hjust=1))
		p=p+theme(legend.position="bottom")
	}
	return (p)
}
#ggsave(p,file='C:/Users/wangzy/Desktop/test.pdf',width=10,height=1.5)













#========================================================================================================================
# 2020-12-9
# another way to plot inferCNV 
#========================================================================================================================

CalCulateCNV=function(infercnv_obj,hpdat,SampleName,chrinfo,mean_hpdat){
 options(stringsAsFactors=F)
 library(infercnv,lib.loc='/public/workspace/lily/R/x86_64-pc-linux-gnu-library/3.6.0/')
 library(TxDb.Hsapiens.UCSC.hg19.knownGene)
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

# plotCNV=function(SEG,limit=c(-0.1,0.1)){
#  library(ggplot2)
#  idx = rev(table(unique(SEG[,c(1,2)])$chr))
#  idx[idx>=2]=2
#  height=c(0)
#  for(i in 1:length(idx)){height=c(height,sum(idx[1:i]))}
#  bk=c()
#  for(i in 1:length(idx)){bk=c(bk,sum(idx[1:i])-idx[i])}
#  SEG = arrange(SEG,Sample,desc(chr),desc(pos))
#  MIN=(limit[1]-min(SEG$seg))/(max(SEG$seg)-min(SEG$seg))
#  MAX=(limit[2]-min(SEG$seg))/(max(SEG$seg)-min(SEG$seg))
 
#  p=ggplot() + 
#   geom_bar(data = SEG, aes(x = Sample,y=1,fill = seg,group=chr),stat = "identity", lwd = 0,color='NA')
#   p=p+geom_hline(aes(yintercept=height), colour="black",lwd=0.5)
#   p=p+coord_flip()+scale_fill_gradientn(colours = c('steelblue',"white","white",'darkred'),
#   values=c(0,MIN,MAX,1))
#   p=p+scale_y_continuous(name = "chromosome",breaks = bk,label=rev(levels(SEG$chr)),expand=c(0,0))+
#   scale_x_discrete(name="",expand=c(0,0))
#   p=p+ theme_set(theme_bw())+theme(panel.grid.major=element_line(colour=NA))
#   p=p+theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),panel.grid.major.y = element_blank(),
#    panel.grid.minor.y = element_blank())
#   p=p+theme(axis.text.x = element_text(angle=45,size=7,hjust=1))
#   p=p+theme(legend.position="bottom")
#  return (p)
# }


library(data.table)
library(ggplot2)
library(dplyr)
chrinfo <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
chrinfo$name=gsub('[0-9].*','',chrinfo$name)
ChromosomeList=paste0('chr',1:22)

objfile=Sys.glob('/public/workspace/lily/Lung2Brain/inte7/infercnv/*/run.final.infercnv_obj')
obvfile=Sys.glob('/public/workspace/lily/Lung2Brain/inte7/infercnv/*/infercnv.observations.txt')

#=========================================================================================================
SEG=data.frame()
for(i in 1:length(objfile)){
  infercnv_obj=readRDS(objfile[i]);
  hpdat=read.table(obvfile[i],header=T,row.names=1,sep=' ',stringsAsFactors=F)
  hpdat=signif(hpdat,3)
  #mean_hpdat=apply(hpdat,1,mean)
  if(i==1){
	  mean_hpdat=apply(hpdat,1,function(x){sum(ifelse(x<0.85,-1,ifelse(x>1.15,1,0)))})
  }
  mean_hpdat=apply(hpdat,1,function(x){sum(ifelse(x<0.95,-1,ifelse(x>1.05,1,0)))})
  mean_hpdat=mean_hpdat/ncol(hpdat)
  SampleName=gsub("/public/workspace/lily/Lung2Brain/inte7/infercnv/|/infercnv.observations.txt","",obvfile[i])
  SEG=rbind(SEG,CalCulateCNV(infercnv_obj, hpdat,SampleName,chrinfo,mean_hpdat))
}
targetinfo=data.frame(chr=rep(rep(ChromosomeList,rep(2,length(ChromosomeList))),length(unique(SEG$Sample))),
pos=rep(c('p','q'),length(ChromosomeList)*length(unique(SEG$Sample))),
Sample=rep(unique(SEG$Sample),rep(2*length(ChromosomeList),length(unique(SEG$Sample)))))
SEG=dplyr::left_join(targetinfo,SEG,by=c('chr','pos','Sample'))
SEG[is.na(SEG$seg),'seg']=0
SEG$chr=factor(SEG$chr,levels=rev(paste0('chr',1:22)))
SEG$pos=factor(SEG$pos,levels=c('q','p'))
SEG$Sample <- factor(SEG$Sample,levels=c("A20190305","A20190312","T-Bsc1","BT1296","BT1297","scrBT1431m","scrBT1432m"))





plotCNV=function(SEG,limit=c(-0.3,0.3)){
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

ggsave(p,file='/public/workspace/lily/Lung2Brain/inte7/Fig/infercnv_zy.pdf',useDingbats = F)







#=============================================================================================================================
# 




















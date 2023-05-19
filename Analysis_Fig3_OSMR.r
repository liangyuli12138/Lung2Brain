
# 2022-6-11 
# This program is used to analysis about OSMR
# 1. find OSMR
# use TCGA data
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
gene <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/chr5p_gene.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/TCGA_LUAD_mod.RDS")

# filter some gene
gene.f <- gene[which(gene%in%rownames(dat))]
tmp.res <- data.frame(t(sapply(gene.f,function(x){
    c(cor.test(as.numeric(dat[x,]),mod[,2],method="spearman")$estimate,
    cor.test(as.numeric(dat[x,]),mod[,2],method="spearman")$p.value)
    })))

colnames(tmp.res) <- c("Rho","pvalue")
tmp.res$exp <- unname(sapply(gene.f,function(x){median(as.numeric(dat[x,]))}))
tmp.res <- tmp.res[order(tmp.res$Rho,decreasing=T),]


# filter some value add some cols
tmp.res.f <- tmp.res[which(tmp.res$pvalue<0.05),]
tmp.res.f$group <- "negative"
tmp.res.f$group[which(tmp.res.f$Rho>0)] <- "positive"

# plot result 
library(ggplot2)
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/TCGA_OSMR_plot.pdf",useDingbats=F)
ggplot(tmp.res.f,aes(x=Rho,y=exp,size= -log2(pvalue),color=group))+geom_point()+scale_color_manual(values=c("#2e9df7","#ec2c22"))+
    theme_bw()
dev.off()






# 1.1 use another dataset to do correlation 
# GSE31210
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_V6","intraInvasion_cp"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE31210_OSMR_BMSCor.pdf",useDingbats=F)
plot(mod$BMS_V6_norm,as.numeric(mod$OSMR),main=" GSE31210 BMS OSMR ")
abline(lm(as.numeric(mod$OSMR)~mod$BMS_V6_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$BMS_V6_norm,as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod$BMS_V6_norm,as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()



























# 2. check OSMR expression percentage in scRNA-seq data
#========================================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
dat$OSMR <- ifelse(dat[["RNA"]]@data["OSMR",]>0,"Y","N")
dat$OSM <- ifelse(dat[["RNA"]]@data["OSM",]>0,"Y","N")



library(ggpubr)
#OSMR
tmp.dat <- apply(table(dat$OSMR,dat$celltype.refine),1,function(x){x/sum(x)})
tmp.dat <- data.frame(tmp.dat)
tmp.dat$celltype <- rownames(tmp.dat)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/scRNA_OSMR_celltype_pie.pdf",useDingbats=F)
ggpie(tmp.dat, 'Y',  
    fill = 'celltype',
    palette = c('#377EB8','#910241','#fc636b','#984EA3','#F29403','#aea400','#B2DF8A',"#E41A1C"), 
    label = paste0(round(tmp.dat$Y * 100,1),"%"), 
    lab.pos = 'in', lab.font = c(4, 'white') #设置标签，标签的位置在图的内部，标签的大小为4， 颜色为白色.
) 
dev.off()

# plot in different group 

tumor <- subset(dat,cells=which(dat$celltype.refine=="Tumor"))
tmp.res <- data.frame(t(apply(table(tumor$OSMR,tumor$type_group),2,function(x){x/sum(x)})))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/scRNA_OSMR_type_group_bar.pdf",useDingbats=F)
barplot(tmp.res[,2],names=rownames(tmp.res),ylim=c(0,0.8))
# 0.6738779 0.5600000 0.4328687
dev.off()





# OSM
tmp.dat <- apply(table(dat$OSM,dat$celltype.refine),1,function(x){x/sum(x)})
tmp.dat <- data.frame(tmp.dat)
tmp.dat$celltype <- rownames(tmp.dat)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/scRNA_OSM_celltype_pie.pdf",useDingbats=F)
ggpie(tmp.dat, 'Y',  
    fill = 'celltype',
    palette = c('#377EB8','#910241','#fc636b','#984EA3','#F29403','#aea400','#B2DF8A',"#E41A1C"), 
    label = paste0(round(tmp.dat$Y * 100,1),"%"), 
    lab.pos = 'in', lab.font = c(4, 'white') #设置标签，标签的位置在图的内部，标签的大小为4， 颜色为白色.
) 
dev.off()



myeloid <- subset(dat,cells=which(dat$celltype.refine=="Myeloid"))
tmp.res <- data.frame(t(apply(table(myeloid$OSM,myeloid$type_group),2,function(x){x/sum(x)})))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/scRNA_OSM_type_group_bar.pdf",useDingbats=F)
barplot(tmp.res[,2],names=rownames(tmp.res),ylim=c(0,0.5))
#  0.2348314 0.3089726 0.2594237

dev.off()












# 3. check OMSR with intravasition and EMT signature 
# in TCGA and GSE200563
# TCGA 

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TCGA_LUAD/TCGA_LUAD_exp.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])
cor.test(mod$OSMR,mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,method="spearman")
cor.test(mod$OSM,mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,method="spearman")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/TCGA_OSMR_EMTCor.pdf",useDingbats=F)
plot(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),main=" TCGA EMT ")
abline(lm(as.numeric(mod$OSMR)~mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()

# calculate extravasation signatures
mod <- mod.analyze2(as.matrix(dat),c("intraInvasion_cp","OSMR_OSM"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])
cor.test(mod$OSMR,mod[,2],method="spearman")
cor.test(mod$OSM,mod[,2],method="spearman")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/TCGA_OSMR_ExtravasationCor.pdf",useDingbats=F)
plot(mod[,2],as.numeric(mod$OSMR),main=" TCGA Extra ")
abline(lm(as.numeric(mod$OSMR)~mod[,2]),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod[,2],as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod[,2],as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()



# GSE200563 data 
# just use MLUNG data 
# 2022-6-13 , test result not O.K.
tmp.dat <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_data.RDS")
tmp.sampleinfo <- readRDS("~/metastasis/data/verify/GSE200563/GSE200563_sampleinfo.RDS")
# just filter MLUNG samples
sampleinfo <- tmp.sampleinfo[which(tmp.sampleinfo$histological%in%c("ADC") & tmp.sampleinfo$group%in%c("MLUNG")),]
dat <- tmp.dat[,rownames(sampleinfo)]

# # calculate signature for EMT
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
# mod <- as.data.frame(mod)
# mod$OSMR <- as.numeric(dat["OSMR",])
# mod$OSM <- as.numeric(dat["OSM",])
# cor.test(mod$OSMR,mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,method="spearman")
# cor.test(mod$OSM,mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,method="spearman")



# GSE68465 data 
# low correlation and 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE68465_OSMR_EMTCor.pdf",useDingbats=F)
plot(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),main=" GSE68465 EMT ")
abline(lm(as.numeric(mod$OSMR)~mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()


# calculate extravasation signatures
mod <- mod.analyze2(as.matrix(dat),c("intraInvasion_cp","BMS_V6","OSMR_OSM"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])
cor.test(mod$OSMR,mod[,2],method="spearman")
cor.test(mod$OSM,mod[,2],method="spearman")


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE68465_OSMR_ExtravasationCor.pdf",useDingbats=F)
plot(mod[,2],as.numeric(mod$OSMR),main=" GSE68465 Extra")
abline(lm(as.numeric(mod$OSMR)~mod$intraInvasion_cp_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod[,2],as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod[,2],as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()








# in GSE31210 data calculate 
# 
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),"/public/workspace/lily/MOD_file/HALLMARK/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE31210_OSMR_EMTCor.pdf",useDingbats=F)
plot(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),main=" GSE31210 EMT ")
abline(lm(as.numeric(mod$OSMR)~mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_norm,as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()



# calculate extravasation signatures
mod <- mod.analyze2(as.matrix(dat),c("intraInvasion_cp"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)
mod$OSMR <- as.numeric(dat["OSMR",])
mod$OSM <- as.numeric(dat["OSM",])


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE31210_OSMR_ExtravasationCor.pdf",useDingbats=F)
plot(mod$intraInvasion_cp_norm,as.numeric(mod$OSMR),main=" GSE31210 Extra")
abline(lm(as.numeric(mod$OSMR)~mod$intraInvasion_cp_norm),col="red")
legend("topright",legend=paste0("rho=",cor.test(mod$intraInvasion_cp_norm,as.numeric(mod$OSMR),method="spearman")$estimate,
    " P =",round(cor.test(mod$intraInvasion_cp_norm,as.numeric(mod$OSMR),method="spearman")$p.value,2))
)
dev.off()







































#========================================================================================================================================
# 4. test in metastasis and non metastasis 
# in GSE126548 not OK
# in GSE123902
library(Seurat)
GSE123902 <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Pyscenic/GSE123902/primary_LUAD.RDS")
GSE123902$group <- "Early"
GSE123902$group[which(GSE123902$sample=="LX675")] <- "Advanced"

GSE123902$OSMR <- ifelse(GSE123902[["RNA"]]@data["OSMR",]>0,"Y","N")
GSE123902$OSM <- ifelse(GSE123902[["RNA"]]@data["OSM",]>0,"Y","N")
apply(table(GSE123902$OSMR,GSE123902$group),2,function(x){x/sum(x)})

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE123902_OSMR_percentage.pdf",useDingbats=F)
barplot(c(0.5786802,0.3443396),names=c("Advanced","Early"),ylab="percentage of OSMR positive tumor cells")
# fisher.test(table(GSE123902$OSMR,GSE123902$group)) P <0.001
dev.off()














# GSE58355 verify 
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/GSE58355/GSE58355_exp.RDS")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(dat),c("BMS_V6"),"/public/workspace/lily/MOD_file/",permN=0)
# mod <- as.data.frame(mod)
sampleinfo <- data.frame(row.names=colnames(dat),sample=colnames(dat),stringsAsFactors=F,
    group=c(rep("Primary",6),rep("CTC",6),rep("Metastasis",6)),
    time=c(rep("day2",3),rep("day25",3),rep("day10",3),rep("day25",3),rep("day10",3),rep("day25",3))
)
sampleinfo$time_group <- paste0(sampleinfo$group,"_",sampleinfo$time)
sampleinfo$OSMR <- as.numeric(dat["OSMR",])

# set factor
sampleinfo$time_group <- factor(sampleinfo$time_group,levels=c("Primary_day2","Primary_day25","CTC_day10","CTC_day25","Metastasis_day10","Metastasis_day25"))
sampleinfo$group <- factor(sampleinfo$group,levels=c("Primary","CTC","Metastasis"))

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE58355_OSMR.pdf",useDingbats=F)
boxplot(OSMR~time_group,data=sampleinfo)
beeswarm::beeswarm(OSMR~time_group,data=sampleinfo,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

boxplot(OSMR~group,data=sampleinfo)
beeswarm::beeswarm(OSMR~group,data=sampleinfo,col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

dev.off()










# 5. survvival verfiy 
# GSE68465
# check result not OK
# library(survminer)
# library(survival)

# dat<- readRDS("/public/workspace/lily/metastasis/data/verify/GSE68465/GSE68465_LUAD_expr.RDS")
# sampleinfo <- readRDS('/public/workspace/lily/metastasis/data/verify/GSE68465/info.RDS')
# sampleinfo <- sampleinfo[colnames(dat),]
# sampleinfo$OS_status[which(sampleinfo$Status=="Dead")] <- 1

# sampleinfo$group <- ifelse(as.numeric(dat["OSMR",])>quantile(as.numeric(dat["OSMR",]),0.67),"High",
#     ifelse(as.numeric(dat["OSMR",])<quantile(as.numeric(dat["OSMR",]),0.33),"Low","Medium"))

# sampleinfo$mths_to_last_clinical_assessment <- as.numeric(as.vector(sampleinfo$mths_to_last_clinical_assessment))
# surv <- Surv(sampleinfo$mths_to_last_clinical_assessment,sampleinfo$OS_status)
# km <- survfit(surv~group,data=sampleinfo)

# pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE68465_OSMR_surv.pdf",useDingbats=F)
# ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',
# 	title="GSE68465",xlab=" Overall survival (months)")+labs(title="OSMR_expression")
# dev.off()









# GSE31210 survival 
dat <- readRDS("~/metastasis/data/verify/GSE31210/GSE31210_expr.RDS")
sampleinfo <- readRDS("~/metastasis/data/verify/GSE31210/sampleinfo.RDS")
sampleinfo <- sampleinfo[colnames(dat),]

sampleinfo$group <- ifelse(as.numeric(dat["OSMR",])>quantile(as.numeric(dat["OSMR",]),0.67),"High",
    ifelse(as.numeric(dat["OSMR",])<quantile(as.numeric(dat["OSMR",]),0.33),"Low","Medium"))

surv <- Surv(as.numeric(sampleinfo$RFS.time),as.numeric(sampleinfo$RFS))
# surv <- Surv(as.numeric(sampleinfo$OS.time),as.numeric(sampleinfo$OS))
km <- survfit(surv~group,data=sampleinfo)


pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GSE31210_OSMR_surv.pdf",useDingbats=F)
ggsurvplot(km,pval=T,palette=c(rgb(251,176,59,maxColorValue = 255),rgb(140,198,63,maxColorValue = 255),'grey'),surv.median.line='hv',
	title="GSE31210",xlab=" RFS (months)")+labs(title="OSMR_expression")
dev.off()












# 2022-8-16
# GSE162198
# OSM in Macrophage 
#=====================================================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE162698/GSE162698_tpm.RDS")

plotdat <- data.frame(
    OSM=c(dat["OSM","Sp1_M0"],dat["OSM","Sp2_M0"],dat["OSM","Sp3_M0"],
        dat["OSM","Sp1_TAMlike"],dat["OSM","Sp2_TAMlike"],dat["OSM","Sp3_TAMlike"]),
    Sample=c("Sp1","Sp2","Sp3","Sp1","Sp2","Sp3"),
    Group=c("M0","M0","M0","TAM","TAM","TAM")
    )

library(ggplot2)
p <- ggplot(plotdat, aes(x = Group, y = OSM)) +
geom_boxplot(aes(fill = Group), show.legend = FALSE, width = 0.6) +  #绘制箱线图
scale_fill_manual(values = c('#FE7280', '#AC88FF')) +  #箱线图的填充色
geom_point(size = 2) +  #绘制样本点
geom_line(aes(group = Sample), color = 'black', lwd = 0.5) +  #绘制配对样本间连线
##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black', size = 1), panel.background = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
    axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 20, color = 'black')) +
labs(x = 'Macrophage', y = 'expression of OSM', title = 'OSM', subtitle = 'M0 vs. TAM')

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GS162198_OSM.pdf",useDingbats=F)
p
dev.off()






# GSE116946
# OSM in Macrophage 
#=====================================================================================================================================================

dat<- readRDS("~/metastasis/data/verify/GSE116946/GSE116946_dat.RDS")
sampleinfo <- readRDS("~/metastasis/data/verify/GSE116946/GSE116946_sampleinfo.RDS")

sampleinfo$OSM <- as.numeric(dat["OSM",])
sampleinfo.f <- sampleinfo[which(sampleinfo$SID %in% names(which(table(sampleinfo$SID)>1))),]

library(ggplot2)
p <- ggplot(sampleinfo, aes(x = Group, y = OSM)) +
geom_boxplot(aes(fill = Group), show.legend = FALSE, width = 0.6) +  #绘制箱线图
scale_fill_manual(values = c('#FE7280', '#AC88FF')) +  #箱线图的填充色
geom_point(size = 2) +  #绘制样本点
geom_line(aes(group = SID), color = 'black', lwd = 0.5) +  #绘制配对样本间连线
##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black', size = 1), panel.background = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
    axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 20, color = 'black')) +
labs(x = 'Macrophage', y = 'expression of OSM', title = 'OSM', subtitle = 'NTAM vs. TAM')

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GS116946_OSM.pdf",useDingbats=F)
p
dev.off()







# GSE162669
# OSM in Macrophage 
#=====================================================================================================================================================

dat <- readRDS("~/metastasis/data/verify/GSE162669/GSE162669_tpm.RDS")

plotdat <- data.frame(
    row.names=colnames(dat),Name=colnames(dat)
    )
plotdat$Group <- sapply(strsplit(as.character(plotdat$Name),"_"),function(x){x[2]})
plotdat$Sample <- sapply(strsplit(as.character(plotdat$Name),"_"),function(x){paste0(x[1],"_",x[3])})
plotdat$OSM <- as.numeric(dat["OSM",])

library(ggplot2)
p <- ggplot(plotdat, aes(x = Group, y = OSM)) +
geom_boxplot(aes(fill = Group), show.legend = FALSE, width = 0.6) +  #绘制箱线图
scale_fill_manual(values = c('#FE7280', '#AC88FF')) +  #箱线图的填充色
geom_point(size = 2) +  #绘制样本点
geom_line(aes(group = Sample), color = 'black', lwd = 0.5) +  #绘制配对样本间连线
##以下是ggplot2的主题设置，修改边框、背景、标题、字体等
theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black', size = 1), panel.background = element_blank(), 
    plot.title = element_text(size = 20, hjust = 0.5), plot.subtitle = element_text(size = 15, hjust = 0.5), 
    axis.text = element_text(size = 20, color = 'black'), axis.title = element_text(size = 20, color = 'black')) +
labs(x = 'Macrophage', y = 'expression of OSM', title = 'OSM', subtitle = 'M0 vs. TAM')



pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GS162669_OSM.pdf",useDingbats=F)
p
dev.off()










# GSE100412
# OSMR in Tumor cells
#=====================================================================================================================================================
dat <- readRDS("~/metastasis/data/verify/GSE100412/GSE100412_FPKM.RDS")
sampleinfo <- readRDS("~/metastasis/data/verify/GSE100412/GSE100412_sampleinfo.RDS")

sampleinfo$Osm <- as.numeric(dat["Osm",])
sampleinfo$Osmr <- as.numeric(dat["Osmr",])
sampleinfo$Cellline <- gsub("_[0-9]","",sampleinfo$Sample)

wilcox.test(Osmr~Group,data=sampleinfo[which(sampleinfo$Group%in%c("InVitro","InVivo")),])
aggregate(Osmr~Group,data=sampleinfo[which(sampleinfo$Group%in%c("InVitro","InVivo")),],FUN=median)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/GS100412_Osmr.pdf",useDingbats=F)
boxplot(log2(Osmr)~Group,data=sampleinfo[which(sampleinfo$Group%in%c("InVitro","InVivo")),])
beeswarm::beeswarm(log2(Osmr)~Group,data=sampleinfo[which(sampleinfo$Group%in%c("InVitro","InVivo")),],col = 4, pch = 16,main = 'beeswarm + bxplot',add = TRUE)

dev.off()












# 2022-11-03
# try to find infercnv groups [chr5p] correlation with BMS score 
# this run in .3
# 感觉要LCBM/nMLUAD 分开做，不然数据的批次会受到影响
#=====================================================================================================================================================
library(Seurat)
library(infercnv)
# LCBM 

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte_BMS_mod.RDS")
all(colnames(dat)==rownames(mod))
mod$BMS <- "Medium"
mod$BMS[which(mod$BMS_V6_norm>quantile(mod$BMS_V6_norm,0.75))] <- "BMSH"
mod$BMS[which(mod$BMS_V6_norm<quantile(mod$BMS_V6_norm,0.25))] <- "BMSL"
dat$BMS_group <- mod$BMS

sub.dat <- subset(dat,cells=which(dat$BMS_group%in%c("BMSH","BMSL")))
# get T cells for reference
adat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
sub.t <- subset(adat,cells=which(adat$type_group=="LCBM" & adat$celltype=="Tcell"))

sub.dat$CNVtype <- sub.dat$BMS_group
sub.t$CNVtype <- sub.t$celltype
sub <- merge(sub.dat,sub.t)

# run infercnv

DefaultAssay(sub) <- "RNA"

setwd("/public/workspace/lily/Lung2Brain/Version6/InferCNV/LCBM/")
# respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])

write.table(sub$CNVtype,"LCBM_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub@assays$RNA@counts)
write.table(count,"LCBM_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="LCBM_count_exp.txt",
        annotations_file="LCBM_cell_info.txt",
        delim="\t",
        gene_order_file="/public/workspace/lily/REF/INDEX-hg38/hg38_position_pure.txt",
        ref_group_names=c("Tcell"))
            
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="/public/workspace/lily/Lung2Brain/Version6/InferCNV/LCBM/", 
                            cluster_by_groups=T, 
                            plot_steps=F,
                            no_prelim_plot = TRUE,
                            num_threads=20, #big
                            no_plot=F ,
                            output_format = "pdf" # maybe can more quick 
                            # used for final scaling to fit range (0,2) centered at 1.
                            )



# LCBM 所有肿瘤细胞自己做一次
# 2022-11-14
# run in .5
#=====================================================================================================================================================
library(Seurat)
library(infercnv)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/LCBM_Tumor_inte.RDS")
# get T cells for reference
adat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
sub.t <- subset(adat,cells=which(adat$type_group=="LCBM" & adat$celltype=="Tcell"))

dat$CNVtype <- "Tumor"
sub.t$CNVtype <- sub.t$celltype

sub <- merge(dat,sub.t)

# run infercnv
DefaultAssay(sub) <- "RNA"

setwd("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/LCBM/")
# respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])

write.table(sub$CNVtype,"LCBM_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub@assays$RNA@counts)
write.table(count,"LCBM_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="LCBM_count_exp.txt",
        annotations_file="LCBM_cell_info.txt",
        delim="\t",
        gene_order_file="/public/workspace/lily/REF/INDEX-hg38/hg38_position_pure.txt",
        ref_group_names=c("Tcell"))
            
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/LCBM/", 
                            cluster_by_groups=T, 
                            plot_steps=F,
                            no_prelim_plot = TRUE,
                            num_threads=40, #big
                            no_plot=F ,
                            output_format = "pdf" # maybe can more quick 
                            # used for final scaling to fit range (0,2) centered at 1.
                            )







#=====================================================================================================================================================
library(Seurat)
library(infercnv)
# LUAD
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
subdat <- subset(dat,cells=which(dat$type_group=="nMLUAD"))
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(subdat[['RNA']]@data),c("BMS_update"),"/public/workspace/lily/Lung2Brain/Version6/Data/",permN=0)
# mod <- as.data.frame(mod)
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_Tumor_BMS_mod.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/nMLUAD_Tumor_BMS_mod.RDS")

all(colnames(subdat)==rownames(mod))
mod$BMS <- "Medium"
mod$BMS[which(mod$BMS_update_norm>quantile(mod$BMS_update_norm,0.75))] <- "BMSH"
mod$BMS[which(mod$BMS_update_norm<quantile(mod$BMS_update_norm,0.25))] <- "BMSL"
subdat$BMS_group <- mod$BMS

sub.dat <- subset(subdat,cells=which(subdat$BMS_group%in%c("BMSH","BMSL")))

# get T cells for reference
adat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
sub.t <- subset(adat,cells=which(adat$type_group=="nMLUAD" & adat$celltype=="Tcell"))

sub.dat$CNVtype <- sub.dat$BMS_group
sub.t$CNVtype <- sub.t$celltype
sub <- merge(sub.dat,sub.t)

# run infercnv

DefaultAssay(sub) <- "RNA"

setwd("/public/workspace/lily/Lung2Brain/Version6/InferCNV/nMLUAD/")
# respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])

write.table(sub$CNVtype,"nMLUAD_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub@assays$RNA@counts)
write.table(count,"nMLUAD_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="nMLUAD_count_exp.txt",
        annotations_file="nMLUAD_cell_info.txt",
        delim="\t",
        gene_order_file="/public/workspace/lily/REF/INDEX-hg38/hg38_position_pure.txt",
        ref_group_names=c("Tcell"))
            
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="/public/workspace/lily/Lung2Brain/Version6/InferCNV/nMLUAD/", 
                            cluster_by_groups=T, 
                            plot_steps=F,
                            no_prelim_plot = TRUE,
                            num_threads=20, #big
                            no_plot=F ,
                            output_format = "pdf" # maybe can more quick 
                            # used for final scaling to fit range (0,2) centered at 1.
                            )










#=====================================================================================================================================================
library(Seurat)
library(infercnv)
# MLUAD
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
subdat <- subset(dat,cells=which(dat$type_group=="MLUAD"))
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod <- mod.analyze2(as.matrix(subdat[['RNA']]@data),c("BMS_update"),"/public/workspace/lily/Lung2Brain/Version6/Data/",permN=0)
# mod <- as.data.frame(mod)
# saveRDS(mod,file="/public/workspace/lily/Lung2Brain/Version6/Data/MLUAD_Tumor_BMS_mod.RDS")
mod <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/MLUAD_Tumor_BMS_mod.RDS")

all(colnames(subdat)==rownames(mod))
mod$BMS <- "Medium"
mod$BMS[which(mod$BMS_update_norm>quantile(mod$BMS_update_norm,0.75))] <- "BMSH"
mod$BMS[which(mod$BMS_update_norm<quantile(mod$BMS_update_norm,0.25))] <- "BMSL"
subdat$BMS_group <- mod$BMS

sub.dat <- subset(subdat,cells=which(subdat$BMS_group%in%c("BMSH","BMSL")))

# get T cells for reference
adat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")
sub.t <- subset(adat,cells=which(adat$type_group=="MLUAD" & adat$celltype=="Tcell"))

sub.dat$CNVtype <- sub.dat$BMS_group
sub.t$CNVtype <- sub.t$celltype
sub <- merge(sub.dat,sub.t)

# run infercnv

DefaultAssay(sub) <- "RNA"

setwd("/public/workspace/lily/Lung2Brain/Version6/InferCNV/MLUAD/")
# respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV/",samplelist[i])

write.table(sub$CNVtype,"MLUAD_cell_info.txt",sep="\t",col.names=F,quote=F)
count<-as.matrix(sub@assays$RNA@counts)
write.table(count,"MLUAD_count_exp.txt",sep="\t",quote=F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="MLUAD_count_exp.txt",
        annotations_file="MLUAD_cell_info.txt",
        delim="\t",
        gene_order_file="/public/workspace/lily/REF/INDEX-hg38/hg38_position_pure.txt",
        ref_group_names=c("Tcell"))
            
infercnv_obj = infercnv::run(infercnv_obj,
                            cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                            out_dir="/public/workspace/lily/Lung2Brain/Version6/InferCNV/MLUAD/", 
                            cluster_by_groups=T, 
                            plot_steps=F,
                            no_prelim_plot = TRUE,
                            num_threads=20, #big
                            no_plot=F ,
                            output_format = "pdf" # maybe can more quick 
                            # used for final scaling to fit range (0,2) centered at 1.
                            )





#####################################################################################################################################################
# analysis-replot
# 2022-11-5

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

plotCNV=function(SEG,limit=c(-0.1,0.1)){
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


library(data.table)
library(ggplot2)
library(dplyr)
chrinfo <- fread("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
           col.names = c("chrom","chromStart","chromEnd","name","gieStain"))
chrinfo$name=gsub('[0-9].*','',chrinfo$name)
ChromosomeList=paste0('chr',1:22)

infercnv_obj=readRDS("/public/workspace/lily/Lung2Brain/Version6/InferCNV/LCBM/run.final.infercnv_obj");
hpdat=read.table("/public/workspace/lily/Lung2Brain/Version6/InferCNV/LCBM/infercnv.observations.txt",header=T,row.names=1,sep=' ',stringsAsFactors=F)

infercnv_obj=readRDS("/public/workspace/lily/Lung2Brain/Version6/InferCNV/nMLUAD/run.final.infercnv_obj");
hpdat=read.table("/public/workspace/lily/Lung2Brain/Version6/InferCNV/nMLUAD/infercnv.observations.txt",header=T,row.names=1,sep=' ',stringsAsFactors=F)

infercnv_obj=readRDS("/public/workspace/lily/Lung2Brain/Version6/InferCNV/MLUAD/run.final.infercnv_obj");
hpdat=read.table("/public/workspace/lily/Lung2Brain/Version6/InferCNV/MLUAD/infercnv.observations.txt",header=T,row.names=1,sep=' ',stringsAsFactors=F)

# 2022-11-14
# in .5
infercnv_obj=readRDS("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/LCBM/run.final.infercnv_obj");
hpdat=read.table("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/LCBM/infercnv.observations.txt",header=T,row.names=1,sep=' ',stringsAsFactors=F)
hpdat=signif(hpdat,3)
mean_hpdat=apply(hpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
mean_hpdat=mean_hpdat/ncol(hpdat)
SampleName="LCBM"
SEG=CalCulateCNV(infercnv_obj, hpdat,SampleName,chrinfo,mean_hpdat)


hpdat=signif(hpdat,3)
# mean_hpdat=apply(hpdat,1,mean)
# highvalue <- mean(mean(as.numeric(hpdat["ACTB",])),mean(as.numeric(hpdat["GAPDH",]))) + 0.1
# lowvalue <- highvalue - 0.2
# diff BMS high and BMS Low cells
BMSHhpdat <- hpdat[,colnames(infercnv_obj@count.data[,infercnv_obj@observation_grouped_cell_indices$BMSH])]
BMSLhpdat <- hpdat[,colnames(infercnv_obj@count.data[,infercnv_obj@observation_grouped_cell_indices$BMSL])]

BMSHmean_hpdat=apply(BMSHhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
BMSHmean_hpdat=BMSHmean_hpdat/ncol(BMSHhpdat)

BMSLmean_hpdat=apply(BMSLhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
BMSLmean_hpdat=BMSLmean_hpdat/ncol(BMSLhpdat)


BMSHSEG=CalCulateCNV(infercnv_obj, BMSHhpdat,SampleName="BMSH",chrinfo,BMSHmean_hpdat)
BMSLSEG=CalCulateCNV(infercnv_obj, BMSLhpdat,SampleName="BMSL",chrinfo,BMSLmean_hpdat)
SEG <- rbind(BMSHSEG,BMSLSEG)

targetinfo=data.frame(chr=rep(rep(ChromosomeList,rep(2,length(ChromosomeList))),length(unique(SEG$Sample))),
pos=rep(c('p','q'),length(ChromosomeList)*length(unique(SEG$Sample))),
Sample=rep(unique(SEG$Sample),rep(2*length(ChromosomeList),length(unique(SEG$Sample)))))
SEG=dplyr::left_join(targetinfo,SEG,by=c('chr','pos','Sample'))
SEG[is.na(SEG$seg),'seg']=0
SEG$chr=factor(SEG$chr,levels=rev(paste0('chr',1:22)))
SEG$pos=factor(SEG$pos,levels=c('q','p'))


# LCBM
# BMSHmean_hpdat=apply(BMSHhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
# BMSLmean_hpdat=apply(BMSLhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
p=plotCNV(SEG,limit=c(-0.1,0.025))
p
pdf("/public/workspace/lily/Lung2Brain/Version6/InferCNV/LCBM_infercnv_BMS.pdf",useDingbats=F)
p
dev.off()

# 2022-11-14
p=plotCNV(SEG,limit=c(-0.1,0.1))
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig3/LCBM_infercnv_All.pdf",useDingbats=F)
p
dev.off()


# MLUAD
# BMSHmean_hpdat=apply(BMSHhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
# BMSLmean_hpdat=apply(BMSLhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
p=plotCNV(SEG,limit=c(-0.1,0.03))
p
pdf("/public/workspace/lily/Lung2Brain/Version6/InferCNV/MLUAD_infercnv_BMS.pdf",useDingbats=F)
p
dev.off()



# nMLUAD
# BMSHmean_hpdat=apply(BMSHhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
# BMSLmean_hpdat=apply(BMSLhpdat,1,function(x){sum(ifelse(x< 0.93,-1,ifelse(x> 1.05,1,0)))})
p=plotCNV(SEG,limit=c(-0.15,0.05))
p
pdf("/public/workspace/lily/Lung2Brain/Version6/InferCNV/nMLUAD_infercnv_BMS.pdf",useDingbats=F)
p
dev.off()













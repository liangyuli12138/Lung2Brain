
# supplemenatry Figure 1
#======================================================================================================================================
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/inte_S16.RDS")

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Supp_Fig1_landscape_type_group.pdf",useDingbats=F)
DimPlot(dat,group.by="type_group",reduction="tsne",raster=T,cols=c("#f7931e","#006837","#6cbc35"),split.by="type_group")
dev.off()

# maybe not show 
pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Supp_Fig1_landscape_sample.pdf",useDingbats=F)
DimPlot(dat,group.by="orig.ident",reduction="tsne",raster=T)
dev.off()





# Feature plot 
# marker <- c("CD3D","CD3E","CD2","PTPRC",
# 			"IGHM","CD79A","MS4A1",
# 			"CD68","FCGR1A","LYZ",
# 			"CLDN5","CDH5","VWF",
#             "COL1A2","THY1","DCN",
#             "MAG","MOG","CNDP1",
# 			"KRT19","EPCAM","CDH1"
# 			)

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Supp_Fig1_featureplot.pdf",useDingbats=F)
FeaturePlot(dat,features=c("CD3D","CD3E","CD2","PTPRC","IGHM","CD79A","MS4A1","CD68"),raster=T,reduction="tsne",ncol=4)
FeaturePlot(dat,features=c("FCGR1A","LYZ","CLDN5","CDH5","VWF","COL1A2","THY1","DCN"),raster=T,reduction="tsne",ncol=4)
FeaturePlot(dat,features=c("MAG","MOG","CNDP1","KRT19","EPCAM","CDH1"),raster=T,reduction="tsne",ncol=4)
dev.off()







# plot CNV score difference in Tumor and epithelial and T cells 

library(Seurat)
library(infercnv)
samplelist<- gsub("\\.RDS$","",grep("*.RDS$",dir("/public/workspace/lily/Lung2Brain/Version6/Prepare_Data/"),value=T))
samplelist <- samplelist[-which(samplelist=="PLCBM2")]
tmp.list1 <- list()
tmp.list2 <- list()
tmp.t <- list()
for(i in 1:length(samplelist)){
    respath <- paste0("/public/workspace/lily/Lung2Brain/Version6/Prepare/InferCNV_hg38/",samplelist[i],"/")
    dat <- readRDS(paste0(respath,"run.final.infercnv_obj"))
    # calculate all cells CNV.score
    tmp.cnv <- apply(dat@expr.data,2,function(x){sum((x-1)^2)})

    # get reference CNV.score
    # cutoff use 75% of reference CNV score 
    ref.cnv <- quantile(tmp.cnv[dat@reference_grouped_cell_indices$Tcell],0.75)
    tmp.t[[i]] <- tmp.cnv[dat@reference_grouped_cell_indices$Tcell]
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

pdf("/public/workspace/lily/Lung2Brain/Version6/Result/Fig1/Supp_Fig1_CNV_score.pdf",useDingbats=F)
boxplot(unlist(tmp.list1),unlist(tmp.list2),unlist(tmp.t),names=c("Tumor like","Epithelial like","T cells"),
    outline=F,main="CNV.score(sum of squares)")
dev.off()

























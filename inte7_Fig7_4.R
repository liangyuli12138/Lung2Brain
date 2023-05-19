


#===========================================================================================================
# 2021-3-12
# Fig4 
# Cytotrace result 
############################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS") 
tumor <- subset(dat,cells=which(dat$maliganant=="tumor"))
tmp.data <- tumor[["RNA"]]@data
tumor$type.TF <- "other"
tumor$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]==0)] <- "CEBPB"
tumor$type.TF[which(tmp.data["CEBPB",]==0&tmp.data["MYBL2",]>0)] <- "MYBL2"
tumor$type.TF[which(tmp.data["CEBPB",]>0&tmp.data["MYBL2",]>0)] <- "CE cells"
load("/public/workspace/lily/Lung2Brain/inte7/res_Data/cytotrace_result.RData")


pdf("/public/workspace/lily/Lung2Brain/inte7/Fig/stemness_Malignant.pdf",useDingbats=F)
boxplot(cytotrace~type,data=tmp.res,FUN=median,outline=F)

library(ggridges)
library(ggplot2)
ggplot(tmp.res) +
  geom_density_ridges_gradient(
    aes(x = cytotrace, y = type, fill = type,height = ..density..), 
    scale = 1, rel_min_height = 0.01) +
  theme_ridges(font_size = 13, grid = TRUE) 


dev.off()



#============================================================================================================






















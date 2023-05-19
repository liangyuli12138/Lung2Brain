#========================================================================================================================
# maybe is Fig6 
# test drug resisitent 
# 
#========================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/metastasis/data/verify/TKI_multiple_Lung/multiple_LB_tumor.RDS")

# run ssgsea 
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat[["RNA"]]@data),c("BMS_test"),"/public/workspace/lily/Lung2Brain/inte7/",permN=1000)
mod <- data.frame(mod)
mod$type <- dat$analysis
saveRDS(mod,file="/public/workspace/lily/Lung2Brain/inte7/res_Data/TKI_lung_mod.RDS")
#seem like there is no difference in PD and naive and PR 
#========================================================================================================================


































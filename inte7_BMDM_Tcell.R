#!/usr/bin/Rscript
# 2021-4-12
# BMDM in GBM and Lung cancer bain metastasis samples
#=============================================================================================================================
library(Seurat)

dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type.refine=="BMDM"&dat$type_group%in%c("GBM","LCBM")))
saveRDS(sub.dat,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_brain.RDS")

#############################################################################################################################
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/BMDM_brain.RDS")









































































































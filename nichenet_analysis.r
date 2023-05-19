
# 2020-12-8 
# Nichenet work 


#=========================================================================================================================
# do some prepare 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/inte_TME.RDS")
Myeloid <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.dat <- subset(dat,cells=which(dat$type%in%c("malignant","Myeloid")))
sub.dat$type.refine <- sub.dat$type 
sub.dat$type.refine[which(colnames(sub.dat)%in%colnames(Myeloid))] <- Myeloid$type.refine
sub.dat$sample <- sub.dat$orig.ident
sub.dat$sample[grep("RD-20180817-001-SR18271",sub.dat$sample)] <- "lesion1"
sub.dat$sample[grep("RD-20180817-002-SR18271",sub.dat$sample)] <- "lesion2"
sub.dat$sample[grep("T-Bsc1",sub.dat$sample)] <- "T_Bsc1"
# test 1 use LCBM to others 
sub.dat.f <- subset(sub.dat,cells=which(sub.dat$type.refine%in%c("malignant","BMDM")))
sub.dat.f$type_group.refine <- "LCBM"
sub.dat.f$type_group.refine[which(sub.dat.f$type_group%in%c("GBM","LC"))] <- "Other"
saveRDS(sub.dat.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/ALL_Malignant_BMDM.RDS")




#=========================================================================================================================
# 2020-12-8
# use da
#=========================================================================================================================
#=========================================================================================================================

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/ALL_Malignant_BMDM.RDS")

# load some database
#==========================================================================================================================
library(nichenetr)
library(Seurat)
library(tidyverse)

ligand_target_matrix = readRDS("/public/workspace/zhumy/ref/nichenet/ligand_target_matrix.rds")
lr_network <- readRDS("/public/workspace/zhumy/ref/nichenet/lr_network.rds")
weighted_networks = readRDS("/public/workspace/zhumy/ref/nichenet/weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
ligand_tf_matrix <- readRDS("/public/workspace/zhumy/ref/nichenet/ligand_tf_matrix.rds")

dat@active.ident <- as.factor(dat$type.refine)
DefaultAssay(dat) <- "RNA"
nichenet_output = nichenet_seuratobj_aggregate(
    seurat_obj = dat, 
    receiver = "malignant", 
    condition_colname = "type_group.refine", 
    condition_oi = "LCBM", 
    condition_reference = "Other", 
    geneset = "up",                                                                        
    sender = "BMDM", 
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "human",
    filter_top_ligands = FALSE,
    expression_pct = 0.05
)


#saveRDS(nichenet_output, "Macrophage2Malignant_WT_output.RDS")

#====================================================================================================================================
# result 
#===================================
library(ComplexHeatmap)
library(tidyr)
library(RColorBrewer)
ligand_target_df <- nichenet_output$ligand_target_df %>%
	pivot_wider(names_from = target, values_from = weight)
ligand_target_df[is.na(ligand_target_df)] <- 0
ligand_target_df <- as.data.frame(ligand_target_df)
rownames(ligand_target_df) <- ligand_target_df$ligand
ligand_target_df <- ligand_target_df[,-1]

ligand_receptor_df <- nichenet_output$ligand_receptor_df %>%
	pivot_wider(names_from = receptor, values_from = weight)
ligand_receptor_df[is.na(ligand_receptor_df)] <- 0
ligand_receptor_df <- as.data.frame(ligand_receptor_df)
rownames(ligand_receptor_df) <- ligand_receptor_df$ligand
ligand_receptor_df <- ligand_receptor_df[,-1]
ligand_receptor_df <- ligand_receptor_df[rownames(ligand_target_df),]

write.csv(nichenet_output$ligand_target_df, "./result/TF/nichenetr/Macro_Malignant_metastasis/Macrophage2Malignant_MvsnonM_ligand_target.csv")
write.csv(nichenet_output$ligand_receptor_df, "./result/TF/nichenetr/Macro_Malignant_metastasis/Macrophage2Malignant_MvsnonM_ligand_receptor.csv")
write.csv(ligand_target_df, "./result/TF/nichenetr/Macro_Malignant_metastasis/Macrophage2Malignant_MvsnonM_ligand_target2.csv")
write.csv(ligand_receptor_df, "./result/TF/nichenetr/Macro_Malignant_metastasis/Macrophage2Malignant_MvsnonM_ligand_receptor2.csv")

































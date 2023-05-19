
#==========================================================================
# analysis Treg to tumor cellphone Db result 
# 2020-12-17
#==========================================================================

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Treg.Tumor <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Treg.malignant")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Treg.malignant") 
	Treg.Tumor[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Treg.Tumor)){
	tmp.id <- c(tmp.id,as.vector(Treg.Tumor[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



Treg.Tumor.res <- cbind(Treg.Tumor[[1]][which(Treg.Tumor[[1]]$id_cp_interaction%in%co_id),],
							Treg.Tumor[[2]][which(Treg.Tumor[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Treg.Tumor)){
	Treg.Tumor.res <- cbind(Treg.Tumor.res,Treg.Tumor[[i]][which(Treg.Tumor[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
Treg.Tumor.res <- Treg.Tumor.res[,c(1,grep("Treg",colnames(Treg.Tumor.res)))]
# add pair gene information 
Communcation.res <- merge(Treg.Tumor.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_malignant_pvalue.RDS")


#===============================================================================================================================
# means value 

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Treg.Tumor <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Treg.malignant")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Treg.malignant") 
	Treg.Tumor[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Treg.Tumor)){
	tmp.id <- c(tmp.id,as.vector(Treg.Tumor[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



Treg.Tumor.res <- cbind(Treg.Tumor[[1]][which(Treg.Tumor[[1]]$id_cp_interaction%in%co_id),],
							Treg.Tumor[[2]][which(Treg.Tumor[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Treg.Tumor)){
	Treg.Tumor.res <- cbind(Treg.Tumor.res,Treg.Tumor[[i]][which(Treg.Tumor[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
Treg.Tumor.res <- Treg.Tumor.res[,c(1,grep("Treg",colnames(Treg.Tumor.res)))]
# add pair gene information 
Communcation.res <- merge(Treg.Tumor.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_malignant_means.RDS")


#=================================================================================================================================================
# calculate result 
#=====================================
pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_malignant_pvalue.RDS")
means <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_malignant_means.RDS")

means.f <- means[which(rownames(means)%in%rownames(pvalue)),]
which.max(apply(means.f,1,function(x){mean(x)}))







































#===========================================================================================================================
# another calculate 
# Communication between Treg and Both T cell 
#==============================================

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Treg.Tumor <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Treg.Both_T")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Treg.Both_T") 
	Treg.Tumor[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Treg.Tumor)){
	tmp.id <- c(tmp.id,as.vector(Treg.Tumor[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



Treg.Tumor.res <- cbind(Treg.Tumor[[1]][which(Treg.Tumor[[1]]$id_cp_interaction%in%co_id),],
							Treg.Tumor[[2]][which(Treg.Tumor[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Treg.Tumor)){
	Treg.Tumor.res <- cbind(Treg.Tumor.res,Treg.Tumor[[i]][which(Treg.Tumor[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
Treg.Tumor.res <- Treg.Tumor.res[,c(1,grep("Treg",colnames(Treg.Tumor.res)))]
# add pair gene information 
Communcation.res <- merge(Treg.Tumor.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_Both_Tumor_pvalue.RDS")


#===============================================================================================================================
# means value 

sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/")
Treg.Tumor <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	# tmp.mean <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/CellphoneDB_Data/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,c("id_cp_interaction","Treg.Both_T")]
	colnames(tmp.value.f)[2] <- paste0(sample_name[i],".Treg.Both_T") 
	Treg.Tumor[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(Treg.Tumor)){
	tmp.id <- c(tmp.id,as.vector(Treg.Tumor[[i]]$id_cp_interaction))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



Treg.Tumor.res <- cbind(Treg.Tumor[[1]][which(Treg.Tumor[[1]]$id_cp_interaction%in%co_id),],
							Treg.Tumor[[2]][which(Treg.Tumor[[2]]$id_cp_interaction%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(Treg.Tumor)){
	Treg.Tumor.res <- cbind(Treg.Tumor.res,Treg.Tumor[[i]][which(Treg.Tumor[[i]]$id_cp_interaction%in%co_id),])
}


# do some modify
Treg.Tumor.res <- Treg.Tumor.res[,c(1,grep("Treg",colnames(Treg.Tumor.res)))]
# add pair gene information 
Communcation.res <- merge(Treg.Tumor.res,tmp.value[,c("id_cp_interaction","interacting_pair","receptor_a","receptor_b")],by="id_cp_interaction")
Communcation.res <- Communcation.res[-which(Communcation.res$interacting_pair=="CALCA_CALCR"),]
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_Both_Tumor_means.RDS")


#=================================================================================================================================================
# calculate result 
#=====================================
pvalue <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_Both_Tumor_pvalue.RDS")
means <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/TF_CCC_res/Treg_Both_Tumor_means.RDS")

means.f <- means[which(rownames(means)%in%rownames(pvalue)),]
which.max(apply(means.f,1,function(x){mean(x)}))
apply(means.f,1,function(x){mean(x)})[order(apply(means.f,1,function(x){mean(x)}))]




















































########################################################################################################################################
# ###### Thu Jan 14 16:04:26 CST 2021
# CCC in BMDM and Treg 
####### Thu Jan 14 19:26:02 CST 2021  change into Treg to BMDM
#=======================================================================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,grep("Treg_.\\.BMDM|id_cp_interaction",colnames(tmp.value))]
	colnames(tmp.value.f) <- paste0(sample_name[i],colnames(tmp.value.f))
	cm.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]][,1]))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]][,1]%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]][,1]%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]][,1]%in%co_id),])
}


# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
# some pair is same id but not same gene 
# ========================================
pair_ann <- tmp.value[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="A20190305id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/BMDM_Treg.pvalue.f.RDS")



##################################################################################
# means 
#=================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,grep("Treg_.\\.BMDM|id_cp_interaction",colnames(tmp.value))]
	colnames(tmp.value.f) <- paste0(sample_name[i],colnames(tmp.value.f))
	cm.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]][,1]))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]][,1]%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]][,1]%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]][,1]%in%co_id),])
}


# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
# some pair is same id but not same gene 
# ========================================
pair_ann <- tmp.value[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="A20190305id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs 
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Treg_BMDM.means.f.RDS")



#########################################################################################################################################
# check result 
pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Treg_BMDM.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Treg_BMDM.means.f.RDS")
# set some filter 
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]


# analysis by group 
############################################################################
tmp.res <- t(apply(means.res,1,function(x){c(mean(x[c(1,3)]),mean(x[c(2,4,5)]),mean(x[c(2,4,5)])/mean(x[c(1,3)]))}))
tmp.res <- tmp.res[order(tmp.res[,3],decreasing=T),]
































#############################################################################################################################
####### Thu Jan 14 20:10:13 CST 2021
# Treg 2 Ehdothelial
# Endothelial 2 Treg
#============================================================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,grep("Endothelial\\.Treg_|id_cp_interaction",colnames(tmp.value))]
	colnames(tmp.value.f) <- paste0(sample_name[i],colnames(tmp.value.f))
	cm.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]][,1]))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]][,1]%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]][,1]%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]][,1]%in%co_id),])
}


# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
# some pair is same id but not same gene 
# ========================================
pair_ann <- tmp.value[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
# pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="A20190305id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Treg_Endothelial.pvalue.f.RDS")


####################################################################################################################################
# means 
#===================================================================================================================================
options(stringsAsFactors=F)
sample_name <- dir("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC/",sample_name[i],"/res/means.txt"),sep="\t",header=T)
	tmp.value.f <- tmp.value[,grep("Endothelial\\.Treg_|id_cp_interaction",colnames(tmp.value))]
	colnames(tmp.value.f) <- paste0(sample_name[i],colnames(tmp.value.f))
	cm.Treg[[i]] <- tmp.value.f
}

tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]][,1]))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 



cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]][,1]%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]][,1]%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]][,1]%in%co_id),])
}

# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("Treg",colnames(cm.Treg.res)))]
# add pair gene information 
# some pair is same id but not same gene 
# ========================================
pair_ann <- tmp.value[,c("id_cp_interaction","interacting_pair")]
# pair_ann[which(pair_ann$id_cp_interaction%in%names(which(table(pair_ann$id_cp_interaction)>1))),]
# pair_ann <- pair_ann[-c(228,233,242,397),]
Communcation.res <- merge(cm.Treg.res,pair_ann,by.x="A20190305id_cp_interaction",by.y="id_cp_interaction")
rownames(Communcation.res) <- Communcation.res$interacting_pair

#==========================
# filter some pairs 
Com.res <- Communcation.res[,grep("Treg",colnames(Communcation.res))]
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Treg_Endothelial.mean.f.RDS")

###################################################################################################################
# check result  
pvalue.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Endothelial_Treg.pvalue.f.RDS")
means.res <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/Endothelial_CCC_res/Endothelial_Treg.mean.f.RDS")
# set some filter 
means.res <- means.res[which(rownames(means.res)%in%rownames(pvalue.res)),]

# analysis by group 
############################################################################
tmp.res <- t(apply(means.res,1,function(x){c(mean(x[c(1,3)]),mean(x[c(2,4,5)]),mean(x[c(2,4,5)])/mean(x[c(1,3)]))}))
tmp.res <- tmp.res[order(tmp.res[,3],decreasing=T),]



































################################################################################################################################################
# 2021-1-21
# run LCBM Treg with BMDM /MG 
#===============================================================================================================================================
library(Seurat)
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/Myeloid/inte9_Myeloid.RDS")
sub.mye <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM")) # all myeloid in LCBM
sub.mye$celltype.refine <- sub.mye$type.refine
tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/inte8_T_cell.RDS")
sub.tcell <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM")) # all T cell in LCBM 

dat <- merge(sub.mye,sub.tcell)

# make cellphoeDB output
################################################################################################################################################
cellphoneDB_input <- function(mat, clusterInfo, expr_outfile, cellinfo_outfile){      
	# mat: Seurat RDS      
	# clusterInfo: group in mat@meta.data      
	count = as.data.frame(mat@assays$RNA@counts)      
	Gene = rownames(mat@assays$RNA@counts)      
	genes = as.data.frame(Gene)      
	count1 <- cbind(genes, count)      
	write.table(count1, expr_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")      
	info <- data.frame(Cell = colnames(mat),           
	cell_type = mat@meta.data[, clusterInfo])      
	write.table(info, cellinfo_outfile,           
	row.names = FALSE, quote = FALSE, sep = "\t")  
}
# change sample name 
dat$sample <- dat$orig.ident
dat$sample[grep("T-Bsc1",dat$sample)] <- "T_Bsc1"
sample_name <- unique(dat$sample)
for(i in 1:length(sample_name)){
    dir.create(paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/",sample_name[i]))
    outpath <- paste0("/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/",sample_name[i],"/")
    tmp <- subset(dat,cells=which(dat$sample==sample_name[i]))
    cellphoneDB_input(tmp,"celltype.refine",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))   
}

# run in shell 
for i in `ls /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/${i}/res
	RESPATH="/public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/${i}/res/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done



############################################################################
# check result 
# just check BMDM and MG with Treg
# 2021-1-23 check MDM with naive T cell 
###### Fri Jan 22 13:50:57 CST 2021
# ~/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/
options(stringsAsFactors=F)
sample_name <- dir("~/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/")
cm.Treg <- list()
for(i in 1:length(sample_name)){
	tmp.value <- read.table(paste0("~/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_Tcell/",sample_name[i],"/res/pvalues.txt"),sep="\t",header=T)
	tmp.value$pair <- paste0(tmp.value$id_cp_interaction,".",tmp.value$interacting_pair)
	tmp.value.f <- tmp.value[,c(ncol(tmp.value),grep("MG.naive.CD4",colnames(tmp.value)))]
	colnames(tmp.value.f) <- c("pair",paste0(sample_name[i],colnames(tmp.value.f)[grep("MG.naive.CD4",colnames(tmp.value.f))]))
	cm.Treg[[i]] <- tmp.value.f
}
##############################################################################
tmp.id <- c()
for(i in 1:length(cm.Treg)){
	tmp.id <- c(tmp.id,as.vector(cm.Treg[[i]][,1]))
}
co_id <- names(which(table(tmp.id)==3)) # get co_id for intersection pair id 

cm.Treg.res <- cbind(cm.Treg[[1]][which(cm.Treg[[1]][,1]%in%co_id),],
							cm.Treg[[2]][which(cm.Treg[[2]][,1]%in%co_id),])
# get co_id pvalue matrix 
for(i in 3:length(cm.Treg)){
	cm.Treg.res <- cbind(cm.Treg.res,cm.Treg[[i]][which(cm.Treg[[i]][,1]%in%co_id),])
}

# do some modify
cm.Treg.res <- cm.Treg.res[,c(1,grep("CD4",colnames(cm.Treg.res)))]
# add pair gene information
Communcation.res <- cm.Treg.res
rownames(Communcation.res) <- Communcation.res$pair
# filter some pairs which are all pvalue >0.05
Com.res <- Communcation.res[,grep("CD4",colnames(Communcation.res))]
Com.res[Com.res==0] <- min(Com.res[Com.res>0])
Com.res[Com.res>0.05] <- 1
Com.res <- -log2(Com.res)
# choose pairs 
Com.res.f <- Com.res[-which(rowSums(Com.res)==0),]

# Save Result 
saveRDS(Com.res.f,file="~/Lung2Brain/TME/Final_11_28/T_cell/CCC_MDM_MG_T_res/MG_naiveCD4.pvalue.f.RDS")




































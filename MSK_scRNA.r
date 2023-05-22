
# this program is used to analysis MSK GSE123902 data 
# 2021-6-23
#===========================================================================================================================
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/MSK_scRNA/GSE123902.RDS")
# FeaturePlot(inte,features=c("EPCAM","EGFR"),label=T)
# show cluster 6 maybe tumor cell beacuse of EPCAM expression 
tumor <- subset(dat,cells=which(dat$seurat_clusters==6))
saveRDS(tumor,file="/public/workspace/lily/Lung2Brain/Version5/MSK_scRNA/GSE123902_tumor.RDS")
ntumor <- subset(dat,cells=colnames(dat)[-which(dat$seurat_clusters==6)])
saveRDS(ntumor,file="/public/workspace/lily/Lung2Brain/Version5/MSK_scRNA/GSE123902_ntumor.RDS")



#===========================================================================================================================
# Now classify T cell and NK cell 

dat <- readRDS("/public/workspace/lily/Lung2Brain/Version5/MSK_scRNA/GSE123902_ntumor.RDS")
DefaultAssay(dat) <- "RNA"
pdf("/public/workspace/lily/Lung2Brain/Version5/MSK_scRNA/T_cell_vlnplot.pdf",width=10)
DefaultAssay(dat) <- "RNA"
VlnPlot(dat,features=c("MS4A1","CD79A","CD79B"),pt.size=0) # B cell 
VlnPlot(dat,features=c("PTPRC","CD3D","CD3E"),pt.size=0) # CD4 T
VlnPlot(dat,features=c("GZMA","IL2","GZMK","GNLY","GZMB","IFNG"),pt.size=0) # cytotoxic T 
VlnPlot(dat,features=c("FOXP3","IL2RA","TGFB1","IKZF2"),pt.size=0) # Treg 
VlnPlot(dat,features=c("TIGIT","PDCD1","LAG3","HAVCR2"),pt.size=0) # exhausted
VlnPlot(dat,features=c("TCF7","SELL","LEF1","CCR7"),pt.size=0) # naive T
VlnPlot(dat,features=c("IRF4", "CREM", "NR4A2"),pt.size=0) # Th17 
VlnPlot(dat,features=c("MAF", "CXCR5", "CXCL13"),pt.size=0) # Th
dev.off()




























































# cellphone DB
#  run in R 
library(Seurat)
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

outpath="/public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/LCBM/"
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LCBM_tumor_myeloid.RDS")
cellphoneDB_input(tmp,"celltype",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))

outpath="/public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/LUAD/"
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/LUAD_tumor_myeloid.RDS")
cellphoneDB_input(tmp,"celltype",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))

outpath="/public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/GBM/"
tmp <- readRDS("/public/workspace/lily/Lung2Brain/Version5/CSOmap/GBM_tumor_myeloid.RDS")
cellphoneDB_input(tmp,"celltype",paste0(outpath,"expr.txt"),paste0(outpath,"cellinfo.txt"))

# run in shell
for i in `ls /public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/`
do
	bytlib load lib64.pool
	bytlib load sqlite3-snapshot
	bytlib load python-3.6.6
	bytlib load cellphonedb-2.1.1

	mkdir -p /public/workspace/lily/Lung2Brain/Version5/Cellphone/Res/${i}/
	RESPATH="/public/workspace/lily/Lung2Brain/Version5/Cellphone/Res/${i}/"
	# run CellphoneDB
	cellphonedb method statistical_analysis /public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/${i}/cellinfo.txt /public/workspace/lily/Lung2Brain/Version5/Cellphone/Data/${i}/expr.txt --threads 8 \
	--output-path=${RESPATH} --counts-data gene_name

done































































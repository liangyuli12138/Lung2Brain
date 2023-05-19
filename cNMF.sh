

# 2022-6-22 
# analysis LCBM tumor cells use cNMF
# 
# make a count matrix in R
# library(Seurat)
# tmp.dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")
# dat <- subset(tmp.dat,cells=which(tmp.dat$type_group=="LCBM"))

# write.table(t(dat@assays$RNA@counts),file="/public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_tumor.tsv",sep="\t",quote=F,col.names=T,row.names=T)

conda activate cNMF
cnmf --help
#step1
cnmf prepare --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ \
--name LCBM_cNMF -c /public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_tumor.tsv \
-k 3 4 5 6 7 8 9 10  --n-iter 100 --seed 123  --total-workers 10 --numgenes 2000



# step2
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 0 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 1 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 2 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 3 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 4 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 5 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 6 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 7 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 8 --total-workers 10
cnmf factorize --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF --worker-index 9 --total-workers 10



#step3
cnmf combine --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF
# rm /public/workspace/lily/Lung2Brain/Version6/cNMF/LCBM_cNMF/cnmf_tmp/LCBM_cNMF.spectra.k_*.iter_*.df.npz





# step4
cnmf k_selection_plot --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF



#step5
# check result 
cnmf consensus --output-dir /public/workspace/lily/Lung2Brain/Version6/cNMF/ --name LCBM_cNMF \
--components 4 --local-density-threshold 0.2 --show-clustering

























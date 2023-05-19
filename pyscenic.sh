
#==============================================================
# pyscenic 
#==============================================================
# conda env list check 

# conda create -n Pyscenic python=3.6.6
# conda activate Pyscenic
# pip install pyscenic
#=============================================================


#===============================================================================================
# run in shell 
#===============================================================================================
# first you have to generate a csv file 
# run in R 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$type_group=="LCBM"&dat$type=="maliganant"))
write.csv(as.matrix(sub.dat[['RNA']]@data),file="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.csv")

#===============================================================================================
# run in shell use pyscenic 
#===============================================================================================

bytlib load python-3.6.6
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/Scenic/arboreto_with_multiprocessing.py \
    -o /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.adjacencies.csv \
    --method grnboost2 \
    --seed 12345 \
    --num_workers 5 \
    /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.csv \
    /public/workspace/zhumy/ref/SCENIC/hs_hgnc_curated_tfs.txt 

# grnboost2_output="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.adjacencies.csv"
# database_fname1="/public/workspace/zhumy/ref/SCENIC/hg19-500bp-upstream-7species.mc9nr.feather"
# database_fname2="/public/workspace/zhumy/ref/SCENIC/hg19-tss-centered-10kb-7species.mc9nr.feather"
# annotations_fname="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
# expression_mtx_fname="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.csv"
# ctx_output="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant_regulons.tsv"

pyscenic ctx \
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.adjacencies.csv \
/public/workspace/zhumy/ref/SCENIC/hg19-500bp-upstream-7species.mc9nr.feather \
/public/workspace/zhumy/ref/SCENIC/hg19-tss-centered-10kb-7species.mc9nr.feather \
--annotations_fname /public/workspace/lily/Lung2Brain/inte7/Pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.csv \
--output /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant_regulons.tsv \
--num_workers 6 \
--mode "custom_multiprocessing"

#=====================================================================================================

pyscenic aucell \
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant.csv /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant_regulons.tsv \
--output /public/workspace/lily/Lung2Brain/inte7/Pyscenic/LCBM_maligant_auc_mtx.tsv \
--num_workers 5 --seed 12345


pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 5 --seed 12345





#=======================================================================================================
# 2020-12-5 
# run all tumor cell in LCBM and LC samples 
#=======================================================================================================
# first you have to generate a csv file 
# run in R 
library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/inte7/Data/inte7_ann.RDS")
sub.dat <- subset(dat,cells=which(dat$type=="maliganant"))
write.csv(t(as.matrix(sub.dat[['RNA']]@data)),file="/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.csv")

#===============================================================================================
# run in shell use pyscenic 
#===============================================================================================

#!/usr/bin/bash

bytlib load python-3.6.6
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/Scenic/arboreto_with_multiprocessing.py \
    -o /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.adjacencies.csv \
    --method grnboost2 \
    --seed 12345 \
    --num_workers 5 \
    /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.csv \
    /public/workspace/zhumy/ref/SCENIC/hs_hgnc_curated_tfs.txt 



pyscenic ctx \
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.adjacencies.csv \
/public/workspace/zhumy/ref/SCENIC/hg19-500bp-upstream-7species.mc9nr.feather \
/public/workspace/zhumy/ref/SCENIC/hg19-tss-centered-10kb-7species.mc9nr.feather \
--annotations_fname /public/workspace/lily/Lung2Brain/inte7/Pyscenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.csv \
--output /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_regulons.tsv \
--num_workers 6 \
--mode "custom_multiprocessing"



pyscenic aucell \
/public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant.csv /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_regulons.tsv \
--output /public/workspace/lily/Lung2Brain/inte7/Pyscenic/ALL_maligant_auc_mtx.tsv \
--num_workers 5 --seed 12345






















####@ parameters
DATABASE="/public/workspace/zhumy/ref/SCENIC"
OUTDIR="/public/workspace/zhumy/CRC2Liver/result/fibroblast/pyscenic"
# grnboost2
expression_mtx_fname="${OUTDIR}/fibroblast.loom"
tfs_fname=${DATABASE}/hs_hgnc_curated_tfs.txt
grnboost2_output="${OUTDIR}/fibroblast.adjacencies.tsv"
# ctx
database_fname1=${DATABASE}/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
database_fname2=${DATABASE}/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather
annotations_fname=${DATABASE}/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
ctx_output="${OUTDIR}/fibroblast.regulons.tsv"
# aucell
aucell_output="${OUTDIR}/fibroblast.auc_mtx.csv"

python ${DATABASE}/arboreto_with_multiprocessing.py \
    $expression_mtx_fname $tfs_fname --method grnboost2 \
    --output $grnboost2_output \
    --num_workers 5 --seed 12345

pyscenic ctx \
    $grnboost2_output \
    $database_fname1 \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $expression_mtx_fname \
    --output $ctx_output \
    --num_workers 5 --mode "custom_multiprocessing" 

pyscenic aucell \
    $expression_mtx_fname $ctx_output \
    --output $aucell_output \
    --num_workers 5 --seed 12345


library(philentropy)
aucmtx <- read_csv("/public/workspace/zhumy/CRC2Liver/result/fibroblast/pyscenic/fibroblast.auc_mtx.csv")
aucmtx1 <- aucmtx[,-1]

rssMat <- sapply(colnames(aucmtx1), function(x) {
    sapply(unique(as.character(fibro@active.ident)), function(y) {
        tmp = rbind(
            value = as.numeric(aucmtx1[[x]]), 
            groups = as.numeric(as.character(ifelse(as.character(fibro@active.ident) == y, 1, 0)))
        )
        1 - JSD(tmp, unit = 'log2', est.prob = "empirical")
    })
})




































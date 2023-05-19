


# run PySCENIC 
# 2022-4-9
# run by APP

library(Seurat)
dat <- readRDS("/public/workspace/lily/Lung2Brain/Version6/Data/Epithelial.RDS")

tumor <- subset(dat,cells=which(dat$celltype.refine=="Tumor"))
DefaultAssay(tumor) <- "RNA"
saveRDS(tumor,file="/public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS")

options(digits=3)
write.csv(t(as.matrix(tumor[['RNA']]@counts)),file="/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/inteS16_tumor.csv")




# shell
#======================================================================================================================================
# Run pySCENIC by pipeline 
# /.bio-apps/pipelines/pySCENIC_bf7ab7f7-ed4a-4157-8a56-f53c41ed8c59/run.sh -j inteS16_tumor \
# -a /public/workspace/lily/REF/scenic_ref/human/hs_hgnc_curated_tfs.txt \
# -b /public/workspace/lily/REF/scenic_ref/human/hg38.recommend.feather.list.txt \
# -c /public/workspace/lily/REF/scenic_ref/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
# -i /public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS \
# -o /public/workspace/lily/Lung2Brain/Version6/PySCENIC/ -t 10


conda activate pyscenic
preprocessing_fname='/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/inteS16_tumor.csv'
tfs_fname='/public/workspace/lily/REF/scenic_ref/human/hs_hgnc_curated_tfs.txt'
grnboost2_output='/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/step1.adjacencies.tsv'

S=''
A=(`cat /public/workspace/lily/REF/scenic_ref/human/hg38.recommend.feather.list.txt`)
for AA in ${A[@]} ; do
  S=$S' '$AA
done

database_fname=$S
annotations_fname='/public/workspace/lily/REF/scenic_ref/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl'
ctx_output='/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/step2.regulons.tsv'
aucell_output='/public/workspace/lily/Lung2Brain/Version6/PySCENIC/tmp/step3.auc_mtx.csv'
num_threads=10
pyscenic grn $preprocessing_fname $tfs_fname -o $grnboost2_output --num_workers $num_threads


pyscenic ctx \
    $grnboost2_output \
    $database_fname \
    --annotations_fname $annotations_fname \
    --expression_mtx_fname $preprocessing_fname \
    --output $ctx_output \
    --num_workers $num_threads --mode "custom_multiprocessing"



pyscenic aucell \
    $preprocessing_fname $ctx_output \
    --output $aucell_output \
    --num_workers $num_threads --seed 12345

conda deactivate



















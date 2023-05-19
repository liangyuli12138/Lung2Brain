#!/usr/bin/sh



# Run pySCENIC by pipeline 
/.bio-apps/pipelines/pySCENIC_bf7ab7f7-ed4a-4157-8a56-f53c41ed8c59/run.sh -j inteS16_tumor \
-a /public/workspace/lily/REF/scenic_ref/human/hs_hgnc_curated_tfs.txt \
-b /public/workspace/lily/REF/scenic_ref/human/hg38.recommend.feather.list.txt \
-c /public/workspace/lily/REF/scenic_ref/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
-i /public/workspace/lily/Lung2Brain/Version6/Data/Tumor.RDS \
-o /public/workspace/lily/Lung2Brain/Version6/PySCENIC/ -t 10

























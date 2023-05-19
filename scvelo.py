
# 2022-6-9
# use scvelo to analysis velocyto
from functools import cache
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata = scv.read("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Pair_loom.h5ad", cache=False)
scv.utils.show_proportions(adata) # show splice and unsplice
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# new mode 
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata,mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata)
scv.tl.latent_time(adata)

# save data 
# https://github.com/theislab/scvelo/issues/255
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write(file="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res.h5ad")



# 2022-6-29
# check result 
adata = scv.read("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/scvelo_res.h5ad",cache=False)
scv.pl.velocity_embedding_stream(adata,color="orig.ident")
scv.pl.velocity_embedding(adata, color="orig.ident",arrow_length=3, arrow_size=2, dpi=300)

scv.pl.scatter(adata, color=['orig.ident','velocity_pseudotime'], cmap='gnuplot',components='1,3'








# 2022-6-20
# use scvelo to analysis velocyto
# for all LCBM 
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata = scv.read("/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/LCBM_loom/Merge_tumor_loom.h5ad", cache=False)
scv.utils.show_proportions(adata) # show splice and unsplice
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# new mode 
scv.tl.recover_dynamics(adata,n_jobs=20)
scv.tl.velocity(adata,mode="dynamical")
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata)
scv.tl.latent_time(adata)

# save data 
# https://github.com/theislab/scvelo/issues/255
adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
adata.write_h5ad(filename="/public/workspace/lily/Lung2Brain/Version6/Data/Velocyto/Res/Merge_tumor_scvelo_res.h5ad")






library(sceasy)
library(Seurat)
dat <- readRDS("./result/integration/integrated_scVI_all.RDS")
DefaultAssay(dat) <- "RNA"
dat1 <- subset(dat, cells = colnames(dat)[dat$celltype %in% c("Hepatobiliary")])
dat1 <- DietSeurat(dat1, counts = TRUE, data = TRUE, scale.data = FALSE)
dat1@meta.data[, which(str_detect(colnames(dat1@meta.data), "DF.classifications"))] <- NULL
dat1@meta.data[, which(str_detect(colnames(dat1@meta.data), "pANN_"))] <- NULL
sceasy::convertFormat(dat1, from="seurat", to="anndata", main_layer = 'counts', outFile='./result/Hep_Cho/hep_cho_scVI.h5ad')




import scanpy as sc
import pandas as pd
import seaborn as sns
adata = sc.read_h5ad("/public/workspace/zhumy/test/20210926_human_liver/result/Hep_Cho/hep_cho_bbknn.h5ad")
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
# adata.raw = adata
# adata = adata[:, adata.var.highly_variable]
# sc.pp.regress_out(adata, ['nCount_RNA'])
# sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.external.pp.bbknn(adata, batch_key='orig.ident')
sc.tl.umap(adata)
#sc.pl.umap(adata, color=['orig.ident'])
sc.tl.leiden(adata, resolution = 2)
#sc.pl.umap(adata, color=['leiden'], save = "_choepi_leiden.pdf", show=False, use_raw=False)
#sc.pl.umap(adata, color=['groups'], save = "_choepi_groups.pdf", show=False, use_raw=False)
#sc.pl.umap(adata, color=['leiden', 'CD3D', 'ANXA4', 'CYP3A7', 'LYZ', 'DCN',], save = "markers_new2.pdf", show=False, use_raw=False)

#adata.obs.groups = adata.obs.groups.astype("category")
sc.tl.embedding_density(adata, groupby='groups')
#adata.obs.groups = adata.obs.groups.astype("str")
sc.pl.embedding_density(adata, groupby='groups', save = "_hep_cho_groups_density.pdf", show=False, use_raw=False)
adata = adata.raw.to_adata()
adata.write_h5ad("./result/Hep_Cho/hep_cho_bbknn.h5ad")






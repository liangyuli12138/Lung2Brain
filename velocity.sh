
# Final lly use 202.195.187.3 conda to trun this result 
# however result is not good 
# this program is used to run volocity 
conda create -n velocyto python=3.6

conda activate velocyto
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install pysam
pip install velocyto

# run bam to 
module load  samtools-1.9

# run in 202.195.187.3
conda activate velocity
velocyto run10x -m /public/workspace/lily/REF/hg19_rmsk.gtf \
    -@ 10 \
    /public/workspace/lily/Lung2Brain/before_6_5/Rawdata/A20190312/ \
    ~/ref/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf







# run python 
# 2021-5-12
#====================================================================================================================================
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization



## for A20190305
adata_obj = scv.read("/public/workspace/lily/Lung2Brain/velocity/LCBM_tumor_veolcyto/A05.loom", cache=True)
adata = scv.read('/public/workspace/lily/Lung2Brain/velocity/LCBM_veolcyto/A20190305.loom', cache=True)
adata_merged_a05 = scv.utils.merge(adata_obj, adata)

# run velocity for A20190305
scv.utils.show_proportions(adata_merged_a05) # show splice and unsplice
adata_merged_a05
scv.pp.filter_and_normalize(adata_merged_a05, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merged_a05, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_merged_a05)
scv.tl.velocity_graph(adata_merged_a05)
scv.tl.tsne(adata_merged_a05)
scv.tl.umap(adata_merged_a05)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap')
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'])
scv.pl.velocity_embedding(adata_merged_a05, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150,color=['seurat_clusters'])



# run velocity for A20190305
adata_obj = scv.read("/public/workspace/lily/Lung2Brain/velocity/LCBM_tumor_veolcyto/A12.loom", cache=True)
adata = scv.read('/public/workspace/lily/Lung2Brain/velocity/LCBM_veolcyto/A20190312.loom', cache=True)
adata_merged_a05 = scv.utils.merge(adata_obj, adata)

# run velocity for A20190305
scv.utils.show_proportions(adata_merged_a05) # show splice and unsplice
adata_merged_a05
scv.pp.filter_and_normalize(adata_merged_a05, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merged_a05, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_merged_a05)
scv.tl.velocity_graph(adata_merged_a05)
scv.tl.tsne(adata_merged_a05)
scv.tl.umap(adata_merged_a05)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap')
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'])
# scv.pl.velocity_embedding(adata_merged_a05, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150,color=['seurat_clusters'])






# run velocity for A20190305
adata_obj = scv.read("/public/workspace/lily/Lung2Brain/velocity/LCBM_tumor_veolcyto/TB.loom", cache=True)
adata = scv.read('/public/workspace/lily/Lung2Brain/velocity/LCBM_veolcyto/T_Bsc1.loom', cache=True)
adata_merged_a05 = scv.utils.merge(adata_obj, adata)

# run velocity for A20190305
scv.utils.show_proportions(adata_merged_a05) # show splice and unsplice
adata_merged_a05
scv.pp.filter_and_normalize(adata_merged_a05, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merged_a05, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_merged_a05)
scv.tl.velocity_graph(adata_merged_a05)
scv.tl.tsne(adata_merged_a05)
scv.tl.umap(adata_merged_a05)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color='seurat_clusters')
# scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap')
# scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'])














#################################################################################################################
# 2021-5-25
# run in 202.195.187.3
# run for CRC

conda activate velocity

python
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization



## for A20190305
adata_obj = scv.read("/public/workspace/zhumy/CRC_liver/result/RD001T.loom")
adata = scv.read('/public/workspace/zhumy/CRC_liver/result/cellranger_result1/RD-20180913-001-SR18282/velocyto/RD-20180913-001-SR18282.loom')
adata_merged_a05 = scv.utils.merge(adata_obj, adata)

# run velocity for A20190305
scv.utils.show_proportions(adata_merged_a05) # show splice and unsplice
adata_merged_a05
scv.pp.filter_and_normalize(adata_merged_a05, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merged_a05, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_merged_a05)
scv.tl.velocity_graph(adata_merged_a05)
scv.tl.tsne(adata_merged_a05)
scv.tl.umap(adata_merged_a05)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap')
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'])
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'], show=True, save="RD001T_dynamic_seurat_clusters.pdf", dpi=300)
scv.pl.velocity_embedding(adata_merged_a05, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150,color=['seurat_clusters'])






## for RD002
import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata_obj = scv.read("/public/workspace/zhumy/CRC_liver/result/cellranger_result1/RD-20180913-002-SR18282/velocyto/RD-20180913-002-SR18282.loom")
adata = scv.read('/public/workspace/zhumy/CRC_liver/result/RD002T.loom')
adata_merged_a05 = scv.utils.merge(adata_obj, adata)

# run velocity for A20190305
scv.utils.show_proportions(adata_merged_a05) # show splice and unsplice
adata_merged_a05
scv.pp.filter_and_normalize(adata_merged_a05, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_merged_a05, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata_merged_a05)
scv.tl.velocity_graph(adata_merged_a05)
scv.tl.tsne(adata_merged_a05)
scv.tl.umap(adata_merged_a05)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap')
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'], show=True, save="RD002T_dynamic_seurat_clusters.pdf", dpi=300)
scv.pl.velocity_embedding_stream(adata_merged_a05, basis='umap', color=['seurat_clusters'])
scv.pl.velocity_embedding(adata_merged_a05, basis='umap', arrow_length=2, arrow_size=1.5, dpi=150,color=['seurat_clusters'])




























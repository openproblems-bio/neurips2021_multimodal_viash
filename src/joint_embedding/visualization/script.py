import logging
import anndata as ad
import scanpy as sc
from matplotlib import pyplot as plt

logging.basicConfig(level=logging.INFO)

## VIASH START

# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.

par = {
    'input_prediction': 'resources_test/joint_embedding/test_resource.prediction.h5ad',
    'input_solution': 'resources_test/joint_embedding/test_resource.solution.h5ad',
    'output_umap': 'resources_test/joint_embedding/test_resource.umap.png',
    'output_emb': 'resources_test/joint_embedding/test_resource.emb.png'
}

## VIASH END

logging.info('Reading `h5ad` files...')
adata = ad.read_h5ad(par['input_prediction'])
adata_solution = ad.read_h5ad(par['input_solution'])

logging.info('Annotating...')
adata.obsm['X_emb'] = adata.X
adata.obs['batch'] = adata_solution.obs['batch'][adata.obs_names]
adata.obs['cell_type'] = adata_solution.obs['cell_type'][adata.obs_names]

logging.info('Nearest neighbours...')
sc.pp.neighbors(adata, use_rep='X_emb')

logging.info('UMAP...')
sc.tl.umap(adata)

logging.info('Plotting...')
# UMAP
sc.pl.umap(adata, color=['batch', 'cell_type'], show=False)
plt.savefig(par['output_umap'])

# Embedding
sc.pl.embedding(adata, basis='X_emb', color=['batch', 'cell_type'])
plt.savefig(par['output_emb'])

# adata.obs['Emb_1'] = adata.X[:, 1]
# adata.obs['Emb_2'] = adata.X[:, 2]
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
#
# sc.pl.scatter(adata, x='Emb_1', y='Emb_2', color='batch', show=False, ax=ax1)
# sc.pl.scatter(adata, x='Emb_1', y='Emb_2', color='cell_type', show=False, ax=ax2)

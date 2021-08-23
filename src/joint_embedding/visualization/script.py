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
    'output': 'resources_test/joint_embedding/test_resource.umap.png'
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
sc.pl.umap(adata, color=['batch', 'cell_type'], show=False)
plt.savefig(par['output'])
